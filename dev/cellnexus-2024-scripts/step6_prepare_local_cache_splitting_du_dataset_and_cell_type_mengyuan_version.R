# Step6
# Group samples by dataset_id, cell_type

# This script sets up a robust and scalable data processing pipeline for single-cell RNA sequencing (scRNA-seq) datasets using the targets package in R, which facilitates reproducible and efficient workflows. Specifically, the code orchestrates the ingestion and preprocessing of multiple SingleCellExperiment objects corresponding to different datasets (dataset_id) and targets (target_name). It leverages high-performance computing resources through the crew package, configuring multiple SLURM-based controllers (tier_1 to tier_4) to handle varying computational loads efficiently.
# 
# The pipeline performs several key steps:
#   
#   1.	Data Retrieval: It reads raw SingleCellExperiment objects for each target, ensuring that only successfully loaded data proceeds further.
# 2.	Normalization: Calculates Counts Per Million (CPM) for each cell to normalize gene expression levels across cells and samples.
# 3.	Data Aggregation: Groups the data by dataset_id and tar_group, then combines the SingleCellExperiment objects within each group into a single object, effectively consolidating the data for each dataset.
# 4.	Metadata Integration: Joins additional metadata, such as cell types, by connecting to a DuckDB database and fetching relevant information from a Parquet file. This enriches the single-cell data with essential annotations.
# 5.	Cell Type Segmentation: Splits the combined SingleCellExperiment objects into separate objects based on cell_type, facilitating downstream analyses that are specific to each cell type.
# 6.	Data Saving with Error Handling: Generates unique identifiers for each cell type within a dataset and saves both the raw counts and CPM-normalized data to specified directories. It includes special handling for cases where a cell type has only one cell, duplicating the data to prevent errors during the saving process.
# 
# By integrating targets, crew, and various data manipulation packages (dplyr, tidyverse, SingleCellExperiment), this script ensures that large-scale scRNA-seq data processing is efficient, reproducible, and capable of leveraging parallel computing resources. It is designed to handle edge cases gracefully and provides a clear framework for preprocessing scRNA-seq data, which is essential for subsequent analyses such as clustering, differential expression, and cell type identification.


library(arrow)
library(dplyr)
library(duckdb)

job::job({
  
  get_file_ids = function(cell_annotation ){
    sample_chunk_df = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue::glue("SELECT * FROM read_parquet('{cell_annotation}')"))
      ) |> 
      # Define chunks
      dplyr::count(dataset_id, sample_id, name = "cell_count") |>  # Ensure unique dataset_id and sample_id combinations
      distinct(dataset_id, sample_id, cell_count) |>  # Ensure unique dataset_id and sample_id combinations
      group_by(dataset_id) |> 
      dbplyr::window_order(dataset_id, cell_count, sample_id) |>  # Ensure order. Note: order cell_count only is not enough because it needs a secondary tie-breaker
      mutate(sample_index = row_number()) |>  # Create sequential index within each dataset
      mutate(sample_chunk = (sample_index - 1) %/% 1000 + 1) |>  # Assign chunks (up to 1000 samples per chunk)
      mutate(sample_pseudobulk_chunk = (sample_index - 1) %/% 250 + 1) |> # Max combination of dataset_id, sample_pseudobulk_chunk and file_id_pseudobulk up to 10000
      mutate(cell_chunk = cumsum(cell_count) %/% 100000 + 1) |> # max 20K cells per sample
      ungroup() 
    
    # Test whether cell_chunk and sample_chunk are unique for this sample
    run_chunk_once <- function(column_name, id) {
      sample_chunk_df |> filter(sample_id == id) |> pull(!!column_name)
    }
    
    sample_chunk_results <- replicate(20, run_chunk_once("sample_chunk", "d6e942a09a140ee8bb6f0c3da8defea4___exp7-human-150well."), simplify = FALSE)
    sample_chunk_identical <- all(sapply(sample_chunk_results[-1], function(x) identical(x, sample_chunk_results[[1]])))
    if (!sample_chunk_identical) {
      stop("Inconsistent sample chunk value was generated in multiple runs, this will lead to file id changes")
    }
    
    cell_chunk_results <- replicate(20, run_chunk_once("cell_chunk",  "d6e942a09a140ee8bb6f0c3da8defea4___exp7-human-150well."), simplify = FALSE)
    cell_chunk_identical <- all(sapply(cell_chunk_results[-1], function(x) identical(x, cell_chunk_results[[1]])))
    if (!cell_chunk_identical) {
      stop("Inconsistent cell chunk value was generated in multiple runs, this will lead to file id changes")
    }
    
    
    tbl(
      dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
      sql(glue::glue("SELECT * FROM read_parquet('{cell_annotation}')"))
    ) |> 
      # Cells in cell_annotation could be more than cells in cell_consensus. In order to avoid NA happens in cell_consensus cell_type column
      mutate(cell_type_unified_ensemble = ifelse(cell_type_unified_ensemble |> is.na(),
                                                 "Unknown",
                                                 cell_type_unified_ensemble)) |>
      
      left_join(sample_chunk_df |> select(dataset_id, sample_chunk, sample_pseudobulk_chunk, cell_chunk, sample_id), copy=TRUE) |> 
      
      # Define chunks
      group_by(dataset_id, sample_chunk, cell_chunk, sample_pseudobulk_chunk, cell_type, sample_id) |>
      summarise(cell_count = n(), .groups = "drop") |>
      group_by(dataset_id, sample_chunk, cell_chunk, cell_type) |>
      dbplyr::window_order(desc(cell_count), sample_id) |> # Important!
      mutate(chunk = cumsum(cell_count) %/% 20000 + 1) |> # max 20K cells per sample
      ungroup() |> 
      as_tibble() |> 
      
      # Single cell file ID
      mutate(file_id_cellNexus_single_cell = 
               glue::glue("{dataset_id}___{sample_chunk}___{cell_chunk}___{cell_type}") |> 
               sapply(digest::digest) |> 
               paste0("___", chunk, ".h5ad") 
      ) |> 
      
      # Pseudobulk file id
      mutate(file_id_cellNexus_pseudobulk = 
               glue::glue("{dataset_id}___{sample_pseudobulk_chunk}") |> 
               sapply(digest::digest) |>
               paste0("___", chunk, ".h5ad"))
    
  }
  
  get_file_ids(
    "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_annotation_2024_Jul.parquet" # MODIFY HERE: input cell annotation parquet
  )  |> 
    write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/file_id_cellNexus_single_cell_2024_Jul.parquet") # MODIFY HERE: output file_id parquet
  
  gc()
  
  con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  dir.create("/vast/scratch/users/shen.m/duckdb_tmp", showWarnings = FALSE) # MODIFY HERE: duckdb temp directory
  
  DBI::dbExecute(
    con,
    "SET temp_directory='/vast/scratch/users/shen.m/duckdb_tmp';" # MODIFY HERE: duckdb temp directory (must match dir.create above)
  )
  
  # Create a view for cell_annotation in DuckDB
  # MODIFY HERE: cell_metadata parquet path inside the SQL string below
  dbExecute(con, "
  CREATE VIEW cell_metadata AS
  SELECT 
    CONCAT(cell_, '___', dataset_id) AS cell_,
    * EXCLUDE (cell_, dataset_id_1, X_umap1, X_umap2, sample_placeholder, cell_type)  -- drop original cell_ and dataset_id_1  
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata.parquet')
")
  
  # MODIFY HERE: cell_annotation parquet path inside the SQL string below
  dbExecute(con, "
  CREATE VIEW hpcell_output_metadata AS
  SELECT * EXCLUDE (
    observation_joinid,
    cell_type_ontology_term_id,
    assay,
    donor_id,
    is_primary_data,
    self_reported_ethnicity,
    tissue,
    azimuth,
    blueprint,
    monaco,
    subsets_Mito_sum,
    subsets_Mito_detected,
    ensemble_joinid,
    cell_type_unified,
    data_driven_ensemble,
    observation_originalid
)
  
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_annotation_2024_Jul.parquet')
")

# MODIFY HERE: file_id parquet path inside the SQL string below (should match the write_parquet output above)
dbExecute(con, "
  CREATE VIEW file_id_cellNexus_single_cell AS
  SELECT 
    dataset_id,
    sample_chunk,
    cell_chunk,
    sample_pseudobulk_chunk,
    cell_type,
    sample_id,
    file_id_cellNexus_single_cell,
    file_id_cellNexus_pseudobulk
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/file_id_cellNexus_single_cell_2024_Jul.parquet')
")

# MODIFY HERE: transformation data frame
dbExecute(con, "
  CREATE VIEW sample_distribution_method_tbl AS
  SELECT 
    sample_2 AS sample_id,
    count_upper_bound,
    feature_thresh AS nfeature_expressed_thresh,
    method_to_apply AS inverse_transform
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/sliced_sample_tbl_2024_Jul.parquet')
")

# Perform the left join and save to Parquet
copy_query <- "
  COPY (
     SELECT 
        cell_metadata.cell_ AS cell_id, -- Rename cell_ to cell_id
        COALESCE(hpcell_output_metadata.alive, FALSE) AS alive, -- Set alive column NULL to FALSE
        cell_metadata.* EXCLUDE (cell_),          -- drop cell_ since it's already aliased as cell_id
        hpcell_output_metadata.* EXCLUDE (cell_, dataset_id, sample_id, alive), -- Deduplicate join keys, and aliased column
        file_id_cellNexus_single_cell.* EXCLUDE (sample_id, dataset_id, cell_type), -- Deduplicate join keys
        sample_distribution_method_tbl.* EXCLUDE (sample_id) -- Deduplicate join keys
      FROM cell_metadata

      LEFT JOIN hpcell_output_metadata
        ON hpcell_output_metadata.cell_ = cell_metadata.cell_
        AND hpcell_output_metadata.dataset_id = cell_metadata.dataset_id
    
      LEFT JOIN file_id_cellNexus_single_cell
        ON file_id_cellNexus_single_cell.sample_id = hpcell_output_metadata.sample_id
        AND file_id_cellNexus_single_cell.dataset_id = hpcell_output_metadata.dataset_id
        AND file_id_cellNexus_single_cell.cell_type = hpcell_output_metadata.cell_type
      
      LEFT JOIN sample_distribution_method_tbl
        ON sample_distribution_method_tbl.sample_id = cell_metadata.sample_id
        
      WHERE cell_metadata.dataset_id NOT IN ('99950e99-2758-41d2-b2c9-643edcdf6d82', '9fcb0b73-c734-40a5-be9c-ace7eea401c9') -- (THESE TWO DATASETS DOESNT contain meaningful data - no observation_joinid etc), thus was excluded in the final metadata.
         
  ) TO  '/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_5_0_mengyuan.parquet' -- MODIFY HERE: output merged metadata parquet (v1_2_2)
  (FORMAT PARQUET, COMPRESSION 'gzip');
"

# Execute the final query to write the result to a Parquet file
dbExecute(con, copy_query)

# Disconnect from the database
dbDisconnect(con, shutdown = TRUE)

print("Done.")
})

# We decided to make cell_id lighter without re-run everything in HPCell pipeline. Here to swap cell_id in the metadata
# cell_map is processed in a separate target script step6_supp_dataset_cell_map.R
job::job({
  con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # Create a view for cell_annotation in DuckDB
  # MODIFY HERE: v1_2_2 merged metadata parquet path inside the SQL string below (should match the COPY TO output above)
  dbExecute(con, "
  CREATE VIEW cell_metadata AS
  SELECT *
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_5_0_mengyuan.parquet')
")
  
  # MODIFY HERE: cell_id dictionary parquet path inside the SQL string below
  dbExecute(con, "
  CREATE VIEW cell_map AS
  SELECT *
  FROM read_parquet('/vast/projects//cellxgene_curated/metadata_cellxgene_mengyuan/file_id_cell_id_dict_v1_1_1_Jul_2024.parquet')
")
  
  # Perform the left join and save to Parquet
  copy_query <- "
  COPY (
     SELECT cell_metadata.*,
            cell_map.* EXCLUDE (cell_id, file_id_cellNexus_single_cell)
      FROM cell_metadata
    
      LEFT JOIN cell_map
        ON cell_metadata.cell_id = cell_map.cell_id
        AND cell_metadata.file_id_cellNexus_single_cell = cell_map.file_id_cellNexus_single_cell

  ) TO  '/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_5_1_mengyuan.parquet' -- MODIFY HERE: output final metadata parquet with new cell IDs (v1_3_2)
  (FORMAT PARQUET, COMPRESSION 'gzip');
"
  
  # Execute the final query to write the result to a Parquet file
  dbExecute(con, copy_query)
  
  # Disconnect from the database
  dbDisconnect(con, shutdown = TRUE)
  
  print("Done.")
  
  
})



# MODIFY HERE: final metadata parquet path used for the targets pipeline (should match the COPY TO output above)
cell_metadata = 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_5_1_mengyuan.parquet')")
  )

library(targets)
library(tidyverse)
store_file_cellNexus = "/vast/scratch/users/shen.m/targets_prepare_database_split_datasets_chunked_1_4_1_single_cell" # MODIFY HERE: targets store directory for this pipeline

tar_script({
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(targets)
  library(tarchetypes)
  library(crew)
  library(crew.cluster)
  
  # Helper (optional) to avoid repetition
  new_elastic <- function(name, mem_gb, time_min, workers, crashes_max, cpus_per_task = 2, backup = NULL) {
    crew_controller_slurm(
      name = name,
      workers = workers,
      crashes_max = crashes_max,
      seconds_idle = 30,
      options_cluster = crew_options_slurm(
        memory_gigabytes_required = mem_gb,
        cpus_per_task = cpus_per_task,
        time_minutes = time_min
      ),
      backup = backup
    )
  }
  
  # Small → large, with fallbacks to the next size up
  elastic_160 <- new_elastic("elastic_160", 160, 60 * 24, workers = 8,  crashes_max = 2)
  elastic_120  <- new_elastic("elastic_120",  120,  60 * 4,  workers = 16, crashes_max = 1, cpus_per_task = 8, backup = elastic_160)
  elastic_80  <- new_elastic("elastic_80",   80,  60 * 4,  workers = 24, crashes_max = 1, cpus_per_task = 8, backup = elastic_120)
  elastic_40  <- new_elastic("elastic_40",   40,  60 * 4,  workers = 32, crashes_max = 1, cpus_per_task = 8, backup = elastic_80)
  elastic_20  <- new_elastic("elastic_20",   20,  60 * 4,  workers = 48, crashes_max = 1, cpus_per_task = 8, backup = elastic_40)
  elastic_10   <- new_elastic("elastic_10",   10, 60 * 4,  workers = 150, crashes_max = 6, cpus_per_task = 8, backup = elastic_20)
  
  elastic_5_minimal   <- new_elastic("elastic_5_minimal",     5, 60 * 4,  workers = 300, crashes_max = 6, cpus_per_task = 8, backup = elastic_10)
  
  
  # Group for targets (small → large)
  controllers <- crew_controller_group(
    elastic_10, elastic_20, elastic_40, elastic_80, elastic_120, elastic_160, elastic_5_minimal
  )
  tar_option_set(
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue", 
    cue = tar_cue(mode = "never"),
    format = "qs",
    #debug = "dataset_id_sct_ea377f6e2d0ae2b7",
    workspace_on_error = TRUE,
    controller = controllers, 
    trust_object_timestamps = TRUE,
    resources = tar_resources(
      crew = tar_resources_crew(controller = "elastic_5_minimal")
    ) 
  )
  
  save_anndata = function(dataset_id_sce, cache_directory){
    
    dir.create(cache_directory, showWarnings = FALSE, recursive = TRUE)
    
    .x = dataset_id_sce |> pull(sce) |> _[[1]]
    .y = dataset_id_sce |> pull(file_id_cellNexus_single_cell) |> _[[1]] |> str_remove("\\.h5ad")
    
    .x |> assays() |> names() = "counts"
    
    # Save the experiment data to the specified counts cache directory
    .x |> save_experiment_data(glue("{cache_directory}/{.y}"))
    
    return(TRUE)  # Indicate successful saving
    
    
  }
  
  # Because they have an inconsistent failure. If I start the pipeline again they might work. Strange.
  insistent_save_anndata <- purrr::insistently(save_anndata, rate = purrr::rate_delay(pause = 60, max_times = 3), quiet = FALSE)
  
  save_anndata_cpm = function(dataset_id_sce, cache_directory){
    
    dir.create(cache_directory, showWarnings = FALSE, recursive = TRUE)
    
    # # Parallelise
    dataset_id_sce |> 
      purrr::transpose() |> 
      lapply(
        FUN = function(x) {
          
          .x = x[[2]]
          .y = x[[1]] |> str_remove("\\.h5ad")
          
          # Check if the 'sce' has only one cell (column)
          if(ncol(assay(.x)) == 1) {
            
            # Duplicate the assay to prevent saving errors due to single-column matrices
            my_assay = cbind(assay(.x), assay(.x))
            # Rename the second column to distinguish it
            colnames(my_assay)[2] = paste0("DUMMY", "___", colnames(my_assay)[2])
            
            cd = colData(.x)
            cd = cd |> rbind(cd)
            rownames(cd)[2] = paste0("DUMMY", "___", rownames(cd)[2])
            
            
            
            .x =  SingleCellExperiment(assay = list( my_assay ) |> set_names(names(assays(.x))[1]), colData = cd) 
          } 
          
          
          # # TEMPORARY FOR SOME REASON THE MIN COUNTS IS NOT 0 FOR SOME SAMPLES
          # .x = HPCell:::check_if_assay_minimum_count_is_zero_and_correct_TEMPORARY(.x, assays(.x) |> names() |> _[1], subset_up_to_number_of_cells = 100)
          
          # CALCULATE CPM
          .x =  SingleCellExperiment(assay = list( cpm = calculateCPM(.x, assay.type = names(assays(.x))[1])), colData = colData(.x)) 
          
          # Save the experiment data to the specified counts cache directory
          .x |> save_experiment_data(glue("{cache_directory}/{.y}"))
          
          return(TRUE)  # Indicate successful saving
        }
        
      )
    
    return("saved")
    
  }
  
  # Because they have an inconsistent failure. If I start the pipeline again they might work. Strange.
  insistent_save_anndata_cpm <- purrr::insistently(save_anndata_cpm, rate = purrr::rate_delay(pause = 60, max_times = 3), quiet = FALSE)
  
  
  # Function to process matrix in vertical slices
  process_matrix_in_slices <- function(h5_matrix, output_filepath, output_filepath_temp, chunk_size = 1000) {
    # Load the HDF5 matrix
    n_rows <- dim(h5_matrix)[1]
    n_cols <- dim(h5_matrix)[2]
    
    if (file.exists(output_filepath)) {
      file.remove(output_filepath)
      cat("Existing output file removed.\n")
    }
    if (file.exists(output_filepath_temp)) {
      file.remove(output_filepath_temp)
      cat("Existing output file removed.\n")
    }
    
    # Create an empty list to hold the slices
    slice_list <- list()
    
    # Loop through the matrix in chunks
    for (start_col in seq(1, n_cols, by = chunk_size)) {
      end_col <- min(start_col + chunk_size - 1, n_cols)
      cat("Processing columns", start_col, "to", end_col, "\n")
      
      # Extract a slice of the matrix
      matrix_slice <- as.matrix(h5_matrix[, start_col:end_col, drop=FALSE])
      
      # Calculate ranks for the slice
      ranked_slice <- singscore::rankGenes(matrix_slice)  %>% `-` (1) 
      
      # Convert the ranked slice to sparse format
      sparse_ranked_slice <- as(ranked_slice, "CsparseMatrix")
      
      # Write the slice to the output HDF5 file
      HDF5Array::writeHDF5Array(
        sparse_ranked_slice,
        filepath = output_filepath_temp,
        name = paste0("rank_", start_col, "_to_", end_col),
        as.sparse = TRUE,
        H5type = "H5T_STD_I32LE"
      ) 
      
      # Store the slice name for later binding
      slice_list[[length(slice_list) + 1]] <- paste0("rank_", start_col, "_to_", end_col)
    }
    
    
    slice_list |> map(~HDF5Array::HDF5Array(output_filepath_temp, name =.x)) |> do.call(cbind, args=_)
    
  }
  
  save_rank_per_cell = function(dataset_id_sce, cache_directory){
    
    dir.create(cache_directory, recursive = TRUE, showWarnings = FALSE)
    
    .x = dataset_id_sce |> pull(sce) |> _[[1]]
    .y = dataset_id_sce |> pull(file_id_cellNexus_single_cell) |> _[[1]] |> str_remove("\\.h5ad")
    
    # Check if the 'sce' has only one cell (column)
    if(ncol(assay(.x)) == 1) {
      
      # Duplicate the assay to prevent saving errors due to single-column matrices
      my_assay = cbind(assay(.x), assay(.x))
      # Rename the second column to distinguish it
      colnames(my_assay)[2] = paste0("DUMMY", "___", colnames(my_assay)[2])
      
      cd = colData(.x)
      cd = cd |> rbind(cd)
      rownames(cd)[2] = paste0("DUMMY", "___", rownames(cd)[2])
      
      
      
      .x =  SingleCellExperiment(assay = list( my_assay ) |> set_names(names(assays(.x))[1]), colData = cd) 
    } 
    
    
    # # TEMPORARY FOR SOME REASON THE MIN COUNTS IS NOT 0 FOR SOME SAMPLES
    # .x = HPCell:::check_if_assay_minimum_count_is_zero_and_correct_TEMPORARY(.x, assays(.x) |> names() |> _[1], subset_up_to_number_of_cells = 100)
    
    print("start ranking")
    
    # CALCULATE rank
    rank_assay = 
      .x |>
      assay() |> 
      
      # This because some datasets are still > 1M cells
      process_matrix_in_slices(
        paste(c(cache_directory, "/", .y, "_rank_matrix.HDF5Array"), collapse = ""), 
        paste(c(cache_directory, "/", .y, "_rank_matrix_temp.HDF5Array"), collapse = ""), 
        chunk_size = 1000
      )
    
    print("creating SCE")
    
    .x =  SingleCellExperiment(assay = list( rank = rank_assay), colData = colData(.x)) 
    
    print("saving")
    
    .x |> save_experiment_data(glue("{cache_directory}/{.y}"))
    
    # Delete the temp file
    file.remove(paste(c(cache_directory, "/", .y, "_rank_matrix_temp.HDF5Array"), collapse = ""))
    
    return(TRUE)  # Indicate successful saving
  }
  
  # Because they have an inconsistent failure. If I start the pipeline again they might work. Strange.
  insistent_save_rank_per_cell <- purrr::insistently(save_rank_per_cell, rate = purrr::rate_delay(pause = 60, max_times = 3), quiet = FALSE)
  
  
  save_anndata_sct = function(dataset_id_sce, cache_directory){
    
    dir.create(cache_directory, showWarnings = FALSE, recursive = TRUE)
    
    if (is.null(dataset_id_sce)) return(NULL)
    
    .x = dataset_id_sce |> pull(sct) |> _[[1]]
    
    # Fix: check is.null BEFORE ncol() to avoid `argument is of length zero`
    if (is.null(.x) || ncol(.x) == 0) return(NULL)
    
    .y = dataset_id_sce |> pull(file_id_cellNexus_single_cell) |> _[[1]] |> str_remove("\\.h5ad")
    
    .x |> assays() |> names() = "sct"
    
    # Wrap save with explicit error logging so the real cause is visible. Strange it shouldnt fail, which passed in debug mode
    tryCatch(
      .x |> save_experiment_data(glue("{cache_directory}/{.y}")),
      error = function(e) {
        message(glue::glue("[save_anndata_sct] FAILED for {.y}: {conditionMessage(e)}"))
        stop(e)
      }
    )
    
    return(TRUE)
  }
  
  # Because they have an inconsistent failure. If I start the pipeline again they might work. Strange.
  insistent_save_anndata_sct <- purrr::insistently(save_anndata_sct, rate = purrr::rate_delay(pause = 60, max_times = 3), quiet = FALSE)
  
  
  cbind_sce_by_dataset_id = function(target_name_grouped_by_dataset_id, file_id_db_file, cell_id_dict, my_store){
    
    my_dataset_id = unique(target_name_grouped_by_dataset_id$dataset_id) 
    my_file_id = unique(target_name_grouped_by_dataset_id$file_id_cellNexus_single_cell) 
    
    file_id_db = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue("SELECT * FROM read_parquet('{file_id_db_file}')"))
      ) |> 
      filter(dataset_id == my_dataset_id) |>
      select(cell_id, sample_id, dataset_id, file_id_cellNexus_single_cell) 
    
    file_id_db = 
      target_name_grouped_by_dataset_id |> 
      left_join(file_id_db, copy = TRUE)
    
    
    dataset_cell_dict = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue("SELECT * FROM read_parquet('{cell_id_dict}')"))
      )  |> 
      filter(file_id_cellNexus_single_cell == my_file_id)
    
    file_id_db = 
      file_id_db |> 
      left_join(dataset_cell_dict, by = c("file_id_cellNexus_single_cell", "cell_id" ), copy=T  )
    
    # Parallelise
    cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))-1
    # Respect R CMD CHECK core limit if set
    if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) {
      cores <- min(cores, 2L)
    }
    # MulticoreParam need to have enough memory to proceed, otherwise reducer error
    bp <- MulticoreParam(workers = cores , progressbar = TRUE)  # Adjust the number of workers as needed
    
    # Begin processing the data pipeline with the initial dataset 'target_name_grouped_by_dataset_id'
    sce_df = 
      file_id_db |> 
      nest(cells = c(cell_id, new_cell_id)) |> 
      # Read raw data for each 'target_name' and store it in a new column 'sce'
      mutate(
        sce = bplapply(
          sce_target_name,
          FUN = function(x) {
            tar_read_raw(x, store = my_store) |> 
              select(.cell, donor_id, dataset_id, sample_id, cell_type) |> 
              mutate(sample_id = as.factor(sample_id)) # lighter
          },  # Read the raw SingleCellExperiment object
          BPPARAM = bp  # Use the defined parallel backend
        )) |> 
      # This should not be needed, but there are some data sets with zero cells 
      filter(!map_lgl(sce, is.null)) |> 
      mutate(sce = map2(sce, cells, ~ {
        
        cell_map <- setNames(.y$new_cell_id, .y$cell_id) 
        
        .x |> filter(.cell %in% names(cell_map)) %>% 
          {
            colnames(.) <- cell_map[colnames(.)]
            .
          }
        
      }, .progress = TRUE))
    
    if(nrow(sce_df) == 0) {
      warning("this chunk has no rows for somereason.")
      return(NULL)
    }
    
    sce_df |> 
      
      # Group the data by 'dataset_id' and 'tar_group' for further summarization
      mutate(sce = map(sce, ~  SingleCellExperiment(assay = assays(.x), colData = colData(.x)) )) |> 
      
      # Combine all 'sce' objects within each group into a single 'sce' object
      group_by(file_id_cellNexus_single_cell) |> 
      summarise( sce =  list(do.call(cbind, args = sce) ),
                 # A step to check missing cells 
                 cells = list(do.call(rbind, args = cells))) 
  }
  
  cbind_sct_by_dataset_id = function(target_name_grouped_by_dataset_id, file_id_db_file, cell_id_dict, my_store){
    
    my_dataset_id = unique(target_name_grouped_by_dataset_id$dataset_id) 
    my_file_id = unique(target_name_grouped_by_dataset_id$file_id_cellNexus_single_cell) 
    
    file_id_db = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue("SELECT * FROM read_parquet('{file_id_db_file}')"))
      ) |> 
      filter(dataset_id == my_dataset_id) |>
      select(cell_id, sample_id, dataset_id, file_id_cellNexus_single_cell) 
    
    file_id_db = 
      target_name_grouped_by_dataset_id |> 
      left_join(file_id_db, copy = TRUE)
    
    
    dataset_cell_dict = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue("SELECT * FROM read_parquet('{cell_id_dict}')"))
      )  |> 
      filter(file_id_cellNexus_single_cell == my_file_id)
    
    file_id_db = 
      file_id_db |> 
      left_join(dataset_cell_dict, by = c("file_id_cellNexus_single_cell", "cell_id"), copy=T  )
    
    # Parallelise
    cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1)) -1
    # Respect R CMD CHECK core limit if set
    if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) {
      cores <- min(cores, 2L)
    }
    # MulticoreParam need to have enough memory to proceed, otherwise reducer error
    bp <- MulticoreParam(workers = cores , progressbar = TRUE)  # Adjust the number of workers as needed  
    
    # Begin processing the data pipeline with the initial dataset 'target_name_grouped_by_dataset_id'
    sct_df = file_id_db |> 
      nest(cells = c(cell_id, new_cell_id)) %>% 
      # Step 1: Read raw data for each 'target_name' and store it in a new column 'sce'
      mutate(
        sct = bplapply(
          sct_target_name,
          FUN = function(x) {
            if (is.na(x)) {
              return(NULL) # because cant get sample_id and dataset_id from NULL sct_matrix
            }
            
            tar_read_raw(x, store = my_store) |>
              select(.cell, donor_id, dataset_id, sample_id, cell_type) |>
              mutate(sample_id = as.factor(sample_id))
            
          },  # Read the raw SingleCellExperiment object
          BPPARAM = bp  # Use the defined parallel backend
        )) |> 
      # This should not be needed, but there are some data sets with zero cells 
      filter(!map_lgl(sct, is.null)) |> 
      mutate(sct = map2(sct, cells, ~ {
        
        cell_map <- setNames(.y$new_cell_id, .y$cell_id) 
        
        .x |> filter(.cell %in% names(cell_map)) %>% 
          {
            colnames(.) <- cell_map[colnames(.)]
            .
          }
        
      }, .progress = TRUE))
    
    if(nrow(sct_df) == 0) {
      warning("this chunk has no rows for somereason.")
      return(NULL)
    }
    
    sct_df |>
      mutate(
        sct = map(sct, \(x) {
          if (is.null(x)) return(NULL)
          SingleCellExperiment(assays = assays(x), colData = colData(x))
        })
      ) |>
      group_by(file_id_cellNexus_single_cell) |>
      summarise(
        sct = {
          scts <- compact(sct)  # drop NULLs inside each group
          
          list(
            if (length(scts) == 0) {
              NULL
            } else {
              
              # A few big samples do not return all features because it reached R limit 2^31-1 in SCTransform
              common_genes <- cellNexus:::check_gene_overlap(scts)
              
              # subset to intersection genes (and keep same order across objects)
              scts2 <- map(scts, \(z) z[common_genes, , drop = FALSE])
              
              do.call(SummarizedExperiment::cbind, scts2)
            }
          )
        },
        cells = list(do.call(rbind, cells)),
        .groups = "drop"
      )
  }
  
  # Because they have an inconsistent failure. If I start the pipeline again they might work. Strange.
  insistent_cbind_sct_by_dataset_id <- purrr::insistently(cbind_sct_by_dataset_id, rate = purrr::rate_delay(pause = 60, max_times = 3), quiet = FALSE)
  
  
  get_dataset_id = function(target_name, my_store){
    # Try reading the target safely (for some failing targets)
    sce = tryCatch(
      tar_read_raw(target_name, store = my_store),
      error = function(e) return(NULL)
    )
    
    # Still need to catch target_name
    if(sce |> is.null()) return(tibble(sample_id = NA_character_, 
                                       dataset_id= NA_character_, 
                                       target_name= !!target_name))
    
    sce |> 
      
      distinct(sample_id, dataset_id) |> mutate(target_name = !!target_name)
  }
  
  create_chunks_for_reading_and_saving = function(dataset_id_sample_id, cell_metadata){
    
    # Solve sample_id mismatches because some end with .h5ad suffix while others dont 
    dataset_id_sample_id |> 
      
      left_join(
        tbl(
          dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
          sql(glue("SELECT * FROM read_parquet('{cell_metadata}')"))
        )   |> 
          distinct(dataset_id, sample_id, sample_chunk, cell_chunk, file_id_cellNexus_single_cell) |> 
          as_tibble(), 
        copy=T
      )
  }
  
  
  cbind_sce_by_dataset_id_get_missing_cells = function(dataset_id_sce){
    
    dataset_id_sce |>
      mutate(
        missing_cells = map2(
          sce, 
          cells, 
          ~{
            cells_in_sce <- .x |> colnames() |> sort()
            
            cells_in_query <- .y$new_cell_id |> unique() |> sort()
            
            # Find differences
            tibble(cell_id = setdiff(cells_in_query, cells_in_sce))
          }
        )
      ) |> 
      select(file_id_cellNexus_single_cell, missing_cells)
    
  }
  
  
  list(
    
    # The input DO NOT DELETE
    tar_target(my_store, "/vast/scratch/users/shen.m/cellNexus/2024-07-01/process_samples_hpcell_target_store", deployment = "main"), # MODIFY HERE: HPCell targets store to read SCEs from
    tar_target(cache_directory, "/vast/scratch/users/shen.m/cellNexus/cellxgene_2024/0.2.1", deployment = "main"), # MODIFY HERE: output cache directory for saved anndata files
    tar_target(
      cell_metadata,
      "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_5_1_mengyuan.parquet", # MODIFY HERE: final metadata parquet (should match the COPY TO output above)
      packages = c( "arrow","dplyr","duckdb")
      
    ),
    
    tar_target(
      cell_id_dict,
      "/vast/projects//cellxgene_curated/metadata_cellxgene_mengyuan/file_id_cell_id_dict_v1_1_1_Jul_2024.parquet", # MODIFY HERE: cell_id dictionary parquet
      packages = c( "arrow","dplyr","duckdb")
    ),
    
    # pre-calculated counts
    tar_target(
      target_name,
      tar_meta(
        starts_with("sce_transformed_"), 
        store = my_store) |> 
        filter(type=="branch") |> 
        pull(name),
      deployment = "main"
    ),
    tar_target(
      dataset_id_sample_id,
      get_dataset_id(target_name, my_store),
      packages = "tidySingleCellExperiment",
      pattern = map(target_name),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_5_minimal")
      )
    ),
    
    # pre-calculated sct
    tar_target(
      sct_target_name,
      tar_meta(
        starts_with("sct_matrix_"), 
        store = my_store) |> 
        filter(type=="branch") |> 
        pull(name),
      deployment = "main"
    ),
    tar_target(
      sct_dataset_id_sample_id,
      get_dataset_id(sct_target_name, my_store),
      packages = "tidySingleCellExperiment",
      pattern = map(sct_target_name),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_5_minimal")
      )
    ),
    
    # join
    tar_target(
      dataset_id_sample_id_target_names,
      dataset_id_sample_id |> left_join(sct_dataset_id_sample_id, by = c("sample_id", "dataset_id"), copy=T) |>
        dplyr::rename(sce_target_name = target_name.x, 
                      sct_target_name = target_name.y),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_5_minimal")
      )
    ),
    
    tar_target(
      target_name_grouped_by_dataset_id,
      create_chunks_for_reading_and_saving(dataset_id_sample_id_target_names, cell_metadata) |> 
        
        # # FOR TESTING PURPOSE ONLY
        # filter(file_id_cellNexus_single_cell %in% c("e8291e33fc98bc21728255f0a5669015___1.h5ad",
        #                                             "4ceb75a970ccdc9aaa9f9e91b931292f___1.h5ad")) |>
        
        group_by(dataset_id, sample_chunk, cell_chunk, file_id_cellNexus_single_cell) |>
        tar_group(),
      iteration = "group",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_5_minimal")
      ), 
      packages = c("arrow", "duckdb", "dplyr", "glue", "targets")
      
    ),
    
    tar_target(
      dataset_id_sce,
      cbind_sce_by_dataset_id(target_name_grouped_by_dataset_id, cell_metadata, cell_id_dict, my_store = my_store),
      pattern = map(target_name_grouped_by_dataset_id),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "scater", "HDF5Array", "arrow", "dplyr", "duckdb",  "BiocParallel", "parallelly"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_40")
      )
    ),
    
    tar_target(
      dataset_id_sct,
      cbind_sct_by_dataset_id(target_name_grouped_by_dataset_id, cell_metadata, cell_id_dict, my_store = my_store),
      pattern = map(target_name_grouped_by_dataset_id),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "scater", "HDF5Array", "arrow", "dplyr", "duckdb",  "BiocParallel", "parallelly"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_120")
      )
    ),
    
    # This target was run for retrieving missing cells analysis only
    tar_target(
      missing_cells_tbl,
      cbind_sce_by_dataset_id_get_missing_cells(dataset_id_sce),
      pattern = map(dataset_id_sce),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "scater", "arrow", "dplyr", "duckdb",  "BiocParallel", "parallelly", "purrr"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_20")
      )
    ),
    
    
    tar_target(
      save_anndata,
      insistent_save_anndata(dataset_id_sce, paste0(cache_directory, "/counts")),
      pattern = map(dataset_id_sce),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_20")
      )
    ),
    
    tar_target(
      saved_dataset_cpm,
      insistent_save_anndata_cpm(dataset_id_sce, paste0(cache_directory, "/cpm")),
      pattern = map(dataset_id_sce),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_20")
      )
    ),
    
    tar_target(
      saved_dataset_rank,
      insistent_save_rank_per_cell(dataset_id_sce, paste0(cache_directory, "/rank")),
      pattern = map(dataset_id_sce),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly", "HDF5Array"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_20")
      )
    ),
    
    tar_target(
      saved_sct,
      save_anndata_sct(dataset_id_sct, paste0(cache_directory, "/sct")),
      pattern = map(dataset_id_sct),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly", "HDF5Array"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_20")
      )
    )
  )
  
}, script = paste0(store_file_cellNexus, "_target_script.R"), ask = FALSE)

job::job({
  
  tar_make(
    script = paste0(store_file_cellNexus, "_target_script.R"), 
    store = store_file_cellNexus, 
    reporter = "summary" #, callr_function = NULL
  )
  
})

missing_cells_tbl = tar_read(missing_cells_tbl, store = store_file_cellNexus) |> 
  unnest(missing_cells)

#missing_cells_tbl |> write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cells_to_remove_in_metadata_Jul_2024.parquet")
missing_cells_tbl <- read_parquet("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cells_to_remove_in_metadata_Jul_2024.parquet")

filtered_cell_metadata = cell_metadata |> anti_join(missing_cells_tbl, by = c("file_id_cellNexus_single_cell",
                                                                              "new_cell_id" = "cell_id"), copy = T)

filtered_cell_metadata |> 
  collect() |> 
  arrow::write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_5_1_filtered_missing_cells_mengyuan.parquet",
                       compression = "zstd") # MODIFY HERE: output parquet after filtering missing cells

