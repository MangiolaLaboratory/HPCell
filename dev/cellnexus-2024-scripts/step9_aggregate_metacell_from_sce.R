library(targets)
# aggregate metacell from metadata and save assays
store_file_cellNexus = "/vast/scratch/users/shen.m/targets_prepare_database_split_datasets_chunked_1_4_0_metacell/" # MODIFY HERE: targets store directory for this pipeline
cell_metadata_path = "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/metadata.2.2.0.parquet" # MODIFY HERE: cell metadata parquet (used to dynamically derive metacell column names)

tar_script({
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(targets)
  library(tarchetypes)
  library(crew)
  library(crew.cluster)
  library(tidySingleCellExperiment)
  library(SingleCellExperiment)
  library(tidyverse)
  library(glue)
  library(digest)
  library(scater)
  library(arrow)
  library(dplyr)
  library(duckdb)
  library(cellNexus)
  library(BiocParallel)
  library(parallelly)
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
    format = "qs",
    
    workspace_on_error = TRUE,
    controller = controllers, 
    trust_object_timestamps = TRUE
    #workspaces = "dataset_id_sce_52dbec3c15f98d66"
  )
  
  get_ids <- function(cell_metadata, metacell_column) {
    metacell_column <- as.character(metacell_column)
    tbl(dbConnect(duckdb::duckdb(), dbdir = ":memory:"), 
        sql(glue("SELECT * FROM read_parquet('{cell_metadata}')"))) |> 
      filter(!is.na(.data[[metacell_column]])) |>
      distinct(file_id_cellNexus_single_cell) |> pull()
  }
  
  get_sce <- function(cell_metadata, id, metacell_column, cache) {
    metacell_column <- as.character(metacell_column)
    tbl(dbConnect(duckdb::duckdb(), 
                  dbdir = ":memory:"), sql(glue("SELECT * FROM read_parquet('{cell_metadata}')"))) |> 
      filter(!is.na(.data[[metacell_column]])) |>
      filter(empty_droplet == FALSE, alive==TRUE, scDblFinder.class!="doublet") |> # because metacell membership ID was pre-calculated after QC
      filter(file_id_cellNexus_single_cell == id) |> 
      select(cell_id, sample_id, dataset_id, donor_id, file_id_cellNexus_single_cell, cell_type, atlas_id, !!metacell_column) |>
      get_single_cell_experiment(cache_directory = cache, repository = NULL) # this assume SCE are not uploaded to cloud
  }
  
  # aggregate_metacell <- function(sce, metacell) {
  #   cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1)) - 1
  #   bp <- MulticoreParam(workers = cores, progressbar = TRUE)
  #   aggregate_metacell <- aggregateAcrossCells(sce, colData(sce)[, c("sample_id", metacell)], BPPARAM = bp)
  #   aggregate_metacell = aggregate_metacell |> mutate(cell_id = paste(sample_id, .data[[metacell]], sep = "___"))
  #   # Assign cell_id to SCE metadata rownames
  #   rownames(colData(aggregate_metacell)) <- aggregate_metacell$cell_id
  #   aggregate_metacell = aggregate_metacell |> select(-contains(".1"))
  #   aggregate_metacell
  #   
  # }
  aggregate_metacell <- function(sce, metacell) {
    
    cores <- max(
      1,
      as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1)) - 1
    )
    
    bp <- BiocParallel::MulticoreParam(
      workers = cores,
      progressbar = TRUE
    )
    
    # aggregate by sample_id + metacell
    agg <- scuttle::aggregateAcrossCells(
      sce,
      ids = colData(sce)[, c("sample_id", metacell), drop = FALSE],
      BPPARAM = bp
    )
    
    # make compact unique cell IDs: <metacell_value>_<index_within_metacell>
    cd <- as.data.frame(SummarizedExperiment::colData(agg))
    
    cd <- cd |>
      dplyr::group_by(.data[[metacell]]) |>
      dplyr::mutate(
        metacell_id = paste0(.data[[metacell]], "_", dplyr::row_number())
      ) |>
      dplyr::ungroup() |> 
      select(-original_cell_)
    
    # put back colData
    SummarizedExperiment::colData(agg) <- S4Vectors::DataFrame(cd)
    
    # use compact .cell as column names / colData rownames
    colnames(agg) <- cd$metacell_id
    rownames(SummarizedExperiment::colData(agg)) <- cd$metacell_id
    
    # optional: remove duplicated helper columns created by aggregation
    keep_cols <- !grepl("\\.1$", colnames(SummarizedExperiment::colData(agg)))
    SummarizedExperiment::colData(agg) <- SummarizedExperiment::colData(agg)[, keep_cols, drop = FALSE]
    
    agg
  }
  
  save_anndata = function(sce, cache_directory) {
    dir.create(cache_directory, showWarnings = FALSE, recursive = TRUE)
    file_id = pull(distinct(sce, file_id_cellNexus_single_cell))
    cellNexus:::save_sce_as_h5ad(sce, glue("{cache_directory}/{file_id}"), mode = "w")
    return(TRUE)
  }
  
  c(
    list(
      tar_target(
        cell_metadata,
        "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/metadata.2.2.0.parquet", # MODIFY HERE: cell metadata parquet (must match cell_metadata_path above)
        deployment = "main",
        packages = c("arrow", "dplyr", "duckdb")
      ),
      tar_target(
        local_cache,
        "/vast/scratch/users/shen.m/cellNexus", # MODIFY HERE: local cache directory containing the single-cell h5ad files (input to get_single_cell_experiment)
        deployment = "main"
      ),
      tar_target(
        save_cache_directory,
        "/vast/scratch/users/shen.m/cellNexus/cellxgene_2024/0.2.0", # MODIFY HERE: output directory where aggregated metacell anndata files are saved
        deployment = "main"
      )
    ),
    tarchetypes::tar_map(
      values = tibble(
        metacell_column = tbl(
          dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
          sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/metadata.2.2.0.parquet')") # MODIFY HERE: cell metadata parquet path inside SQL (must match cell_metadata_path above)
        ) |> select(contains("metacell_")) |> colnames()
      ),
      names = metacell_column,
      unlist = TRUE,
      tar_target(
        file_ids,
        get_ids(cell_metadata, metacell_column) 
        # |>
        #   # TEST PURPOSE ONLY
        #   head(2)
      ),
      tar_target(
        file_id_sce,
        get_sce(cell_metadata, file_ids, metacell_column, local_cache),
        pattern = map(file_ids),
        resources = tar_resources(crew = tar_resources_crew(controller = "elastic_10"))
      ),
      tar_target(
        metacell,
        aggregate_metacell(file_id_sce, metacell_column),
        pattern = map(file_id_sce),
        resources = tar_resources(crew = tar_resources_crew(controller = "elastic_20"))
      ),
      tar_target(
        save_metacell,
        save_anndata(metacell, paste0(save_cache_directory, "/", metacell_column, "/counts")),
        pattern = map(metacell),
        resources = tar_resources(crew = tar_resources_crew(controller = "elastic_20"))
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
#tar_invalidate(names = everything(), store = store_file_cellNexus)
tar_meta(store = store_file_cellNexus) |> filter(!is.na(error)) |> distinct(name, error)
tar_errored(store = store_file_cellNexus)
# With tar_map, target names are suffixed by metacell_column, e.g. file_ids_metacell_4, metacell_metacell_4
# tar_workspace("metacell_metacell_4_<hash>", store = store_file_cellNexus, script = paste0(store_file_cellNexus, "_target_script.R"))
# debugonce(aggregate_metacell)
# aggregate_metacell(file_id_sce, "metacell_4")  # when debugging a specific branch

# Check the number of file id should be created for metacell_2
cache_dir = "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/" # MODIFY HERE: directory used for verification queries below (contains metadata.2.2.0.parquet)
# Define all metacell levels
metacell_levels <- c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536)
metacell_names <- paste0("metacell_", metacell_levels)

# Function to get file IDs that SHOULD be generated (from metadata)
file_ids <- function(metacell_col_name) {
  get_metadata(
    cache_directory = cache_dir,
    cloud_metadata = NULL,
    local_metadata = file.path(cache_dir, "metadata.2.2.0.parquet")
  ) |>
    filter(
      !is.na(.data[[metacell_col_name]]),
      empty_droplet == FALSE,
      alive == TRUE,
      scDblFinder.class != "doublet"
    ) |>
    distinct(file_id_cellNexus_single_cell) |>
    pull(file_id_cellNexus_single_cell) |>
    sort()
}

# Function to get file names that ACTUALLY exist in subfolder
file_names_saved <- function(metacell_name) {
  subfolder <- file.path("~/scratch/cellNexus/cellxgene/01-07-2024/", metacell_name, "counts")
  if (dir.exists(subfolder)) {
    list.files(subfolder, recursive = FALSE) |> sort()
  } else {
    character(0)
  }
}

# Build the summary tibble
result <- tibble(metacell = metacell_names) |>
  mutate(
    ids_to_save  = map(metacell, file_ids,        .progress = "Getting expected file IDs"),
    ids_saved    = map(metacell, file_names_saved, .progress = "Getting actual file names"),
    file_to_save = map_int(ids_to_save, length),
    file_saved   = map_int(ids_saved,   length),
    missing      = map2(ids_to_save, ids_saved, \(expected, actual) setdiff(expected, actual)),
    extra        = map2(ids_to_save, ids_saved, \(expected, actual) setdiff(actual, expected))
  ) |>
  select(metacell, file_to_save, file_saved, missing, extra)

result

# Unit test query lung tissue, metacell 256
lung_metacell_256 = get_metadata(
  cache_directory = cache_dir,
  cloud_metadata = NULL,
  local_metadata = file.path(cache_dir, "metadata.2.2.0.parquet")
) |> 
  filter(!is.na(metacell_256),
         empty_droplet == FALSE,
         alive == TRUE,
         scDblFinder.class != "doublet") |>
  filter(tissue  == "lung") |> 
  get_metacell(cache_directory = cache_dir, cell_aggregation = "metacell_256")

# Cell number in lung_metacell_256 should match the count below:
get_metadata(
  cache_directory = cache_dir,
  cloud_metadata = NULL,
  local_metadata = file.path(cache_dir, "metadata.2.2.0.parquet")
) |> 
  filter(!is.na(metacell_256),
         empty_droplet == FALSE,
         alive == TRUE,
         scDblFinder.class != "doublet") |>
  filter(tissue  == "lung") |> 
  distinct(sample_id, metacell_256, cell_type_unified_ensemble) |> dplyr::count()


