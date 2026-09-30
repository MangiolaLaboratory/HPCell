library(targets)
library(tidyverse)
library(cellNexus)
store_file_cellNexus = "/vast/scratch/users/shen.m/targets_prepare_database_split_datasets_chunked_1_4_1_pseudobulk" # MODIFY HERE: targets store directory for this pipeline
my_store =  "/vast/scratch/users/shen.m/cellNexus/2024-07-01/process_samples_hpcell_target_store" # MODIFY HERE: HPCell targets store to read pseudobulk SCEs from (used throughout)

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
  
  elastic_5_minimal   <- new_elastic("elastic_5_minimal",     5, 60 * 4,  workers = 300, crashes_max = 6, cpus_per_task = 2, backup = elastic_10)
  
  
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
    
    workspace_on_error = TRUE,
    controller = controllers, 
    trust_object_timestamps = TRUE,
    resources = tar_resources(
      crew = tar_resources_crew(controller = "elastic_5_minimal")
    ) 
  )
  
  
  get_dataset_id = function(target_name, my_store){
    sce = tar_read_raw(target_name, store = my_store)
    
    if(sce |> is.null()) return(tibble(sample_id = character(), dataset_id= character(), 
                                       target_name= target_name))
    
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
          distinct(sample_id, sample_pseudobulk_chunk, cell_chunk, 
                   cell_type_unified_ensemble,
                   file_id_cellNexus_pseudobulk) |> 
          as_tibble(), 
        copy=T
      )
  }
  
  
  cbind_sce_by_dataset_id = function(target_name_grouped_by_dataset_id, 
                                     file_id_db_file, my_store){
    
    #my_dataset_id = unique(target_name_grouped_by_dataset_id$dataset_id) 
    my_cell_type = unique(target_name_grouped_by_dataset_id$cell_type_unified_ensemble)
    
    file_id_db = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue("SELECT * FROM read_parquet('{file_id_db_file}')"))
      ) |> 
      dplyr::filter(cell_type_unified_ensemble %in% my_cell_type) |>
      select(sample_id, dataset_id, cell_type_unified_ensemble,
             file_id_cellNexus_pseudobulk) 
    
    
    file_id_db = 
      target_name_grouped_by_dataset_id |> 
      left_join(file_id_db, copy = TRUE)
    
    
    # Parallelise
    cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1)) -1
    # Respect R CMD CHECK core limit if set
    if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) {
      cores <- min(cores, 2L)
    }
    bp <-   MulticoreParam(workers = cores, progressbar = TRUE)
    
    # Begin processing the data pipeline with the initial dataset 'target_name_grouped_by_dataset_id'
    sce_df = 
      file_id_db |> 
      mutate(cell_id = paste(sample_id, cell_type_unified_ensemble, sep = "___")) |>
      nest(cells = cell_id) |> 
      # Step 1: Read raw data for each 'target_name' and store it in a new column 'sce'
      mutate(
        sce = bplapply(
          target_name,
          FUN = function(x) tar_read_raw(x, store = my_store),
          BPPARAM = bp
        )
      ) |>
      
      # This should not be needed, but there are some data sets with zero cells 
      filter(!map_lgl(sce, is.null)) |> 
      
      mutate(sce = map2(sce, cells, ~ .x |>
                          filter(.cell %in% .y$cell_id),
                        
                        .progress = TRUE))
    
    
    
    if(nrow(sce_df) == 0) {
      warning("this chunk has no rows for somereason.")
      return(NULL)
    }
    
    sce_df = sce_df |> 
      
      # THIS SHOULD HAVE BEEN DONE IN THE TRANFORM HPCell
      mutate(sce = map(sce, ~  SingleCellExperiment(assay = assays(.x), colData = colData(.x)) ))
    
    
    # Extra Step 1: Harmonize colData columns - Avoid column name mismatch, force cbind
    all_col_names <- sce_df$sce %>%
      map(~colnames(colData(.x))) %>% 
      unlist() %>% 
      unique()
    
    # Extra Step 2: Standardize colData to have the same columns in each SCE
    sce_df$sce <- map(sce_df$sce, function(sce) {
      current_cols <- colnames(colData(sce))
      missing_cols <- setdiff(all_col_names, current_cols)
      
      if (length(missing_cols) > 0) {
        
        # Fill missing colData columns with NA
        for (col in missing_cols) {
          # Handle sce with empty cells
          if (ncol(sce) == 0)  colData(sce)[, col] <- character(0)
          else if (ncol(sce) > 0) colData(sce)[, col] <- NA
        }
      }
      
      # Ensure the order of columns matches
      colData(sce) <- colData(sce)[, all_col_names]
      return(sce)
    })
    
    sce_df |>
      
      # Step 5: Combine all 'sce' objects within each group into a single 'sce' object
      group_by(file_id_cellNexus_pseudobulk) |>
      summarise( sce =  list(do.call(cbind, args = sce) ) ) 
    
  }
  
  
  
  save_anndata = function(dataset_id_sce, cache_directory){
    
    dir.create(cache_directory, showWarnings = FALSE, recursive = TRUE)
    
    .x = dataset_id_sce |> pull(sce) |> _[[1]]
    .y = dataset_id_sce |> pull(file_id_cellNexus_pseudobulk) |> _[[1]] |> str_remove("\\.h5ad")
    
    .x |> assays() |> names() = "counts"
    
    # Drop list-type columns in colData 
    cd <- colData(.x)
    is_list_col <- vapply(cd, is.list, logical(1))
    colData(.x) <- cd[, !is_list_col, drop = FALSE]
    
    # Check if there is a memory issue 
    assays(.x) <- assays(.x) |> map(DelayedArray::realize)
    
    # Save the experiment data to the specified counts cache directory
    .x |> save_experiment_data(glue("{cache_directory}/{.y}"))
    
    return(TRUE)  # Indicate successful saving
    
  }
  
  # Because they have an inconsistent failure. If I start the pipeline again they might work. Strange.
  insistent_save_anndata <- purrr::insistently(save_anndata, rate = purrr::rate_delay(pause = 60, max_times = 3), quiet = FALSE)
  
  list(
    
    # The input DO NOT DELETE
    tar_target(my_store, "/vast/scratch/users/shen.m/cellNexus/2024-07-01/process_samples_hpcell_target_store", deployment = "main"), # MODIFY HERE: HPCell targets store (must match my_store above)
    tar_target(cache_directory, "/vast/scratch/users/shen.m/cellNexus/cellxgene_2024/0.2.1/pseudobulk", deployment = "main"), # MODIFY HERE: output cache directory for saved pseudobulk anndata files
    tar_target(
      cell_metadata,
      "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_5_1_mengyuan.parquet", # MODIFY HERE: cell metadata parquet (output of step6/step7)
      packages = c( "arrow","dplyr","duckdb")
      
    ),
    tar_target(
      target_name,
      tar_meta(
        starts_with("pseudobulk_se_iterated_"), 
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
        crew = tar_resources_crew(controller = "elastic_10")
      )
    ),
    
    tar_target(
      target_name_grouped_by_dataset_id,
      create_chunks_for_reading_and_saving(dataset_id_sample_id, cell_metadata) |> 
        
        # # FOR TESTING PURPOSE ONLY
        # filter(file_id_cellNexus_pseudobulk %in% c("9722bedfd71d069fe3665b4ae03fbeb9___2.h5ad",
        #                                            "2996bb4263f9fb301d8460f4f0450848___2.h5ad")) |>
        
        group_by(dataset_id,
                 sample_pseudobulk_chunk, 
                 # When using strategy file_id = dataset_id, dont group by cell_chunk as it will result in returning more than one SCEs for the same dataset_id
                 #cell_chunk, 
                 file_id_cellNexus_pseudobulk) |>
        tar_group(),
      iteration = "group",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_5_minimal")
      ), 
      packages = c("arrow", "duckdb", "dplyr", "glue", "targets")
      
    ),
    
    tar_target(
      dataset_id_sce,
      cbind_sce_by_dataset_id(target_name_grouped_by_dataset_id, cell_metadata, my_store = my_store),
      pattern = map(target_name_grouped_by_dataset_id),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "scater", "arrow", "dplyr", "duckdb",  "BiocParallel", "parallelly"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_20")
      )
    ),
    tar_target(
      get_pseudobulk,
      insistent_save_anndata(dataset_id_sce, paste0(cache_directory, "/counts")),
      pattern = map(dataset_id_sce),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly"),
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

