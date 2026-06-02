library(targets)
store = "/vast/scratch/users/shen.m/cellnexus_dataset_cell_map_Nov_2025_v1_0_0_target_store"
tar_script({
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(targets)
  library(tarchetypes)
  library(crew)
  library(crew.cluster)
  
  new_elastic <- function(name, mem_gb, time_min, workers, crashes_max, cpus_per_task = 1, backup = NULL) {
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
  elastic_160 <- new_elastic("elastic_160", 160, 60 * 24, workers = 8,  crashes_max = 2)
  elastic_120  <- new_elastic("elastic_120",  120,  60 * 4,  workers = 16, crashes_max = 1, cpus_per_task = 1, backup = elastic_160)
  elastic_80  <- new_elastic("elastic_80",   80,  60 * 4,  workers = 24, crashes_max = 1, cpus_per_task = 1, backup = elastic_120)
  elastic_40  <- new_elastic("elastic_40",   40,  60 * 4,  workers = 32, crashes_max = 1, cpus_per_task = 1, backup = elastic_80)
  elastic_20  <- new_elastic("elastic_20",   20,  60 * 4,  workers = 48, crashes_max = 1, cpus_per_task = 1, backup = elastic_40)
  elastic_10   <- new_elastic("elastic_10",   10, 60 * 4,  workers = 150, crashes_max = 2, cpus_per_task = 1, backup = elastic_20)
  
  elastic_5_minimal   <- new_elastic("elastic_5_minimal",     5, 60 * 4,  workers = 300, crashes_max = 2, cpus_per_task = 1, backup = elastic_10)
  
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
  
  get_unique_file_ids <- function(cell_metadata){
    tbl(dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue::glue("SELECT * FROM read_parquet('{cell_metadata}')"))) |> 
      distinct(file_id_cellNexus_single_cell) |> pull()
  }
  
  create_file_id_cell_id_dict <- function(cell_metadata, file_id) {
    
    
    tbl(dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue::glue("SELECT * FROM read_parquet('{cell_metadata}')"))) |> 
      filter(file_id_cellNexus_single_cell == file_id) |>
      dbplyr::window_order(cell_id) |> 
      mutate(cell_index = row_number()) |> 
      select(cell_id, file_id_cellNexus_single_cell, 
             new_cell_id = cell_index) |>
      collect()
  }
  
  list(
    tar_target(cell_metadata , "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_metadata_cell_type_consensus_v1_0_0_mengyuan.parquet",
               deployment = "main"),
    tar_target(
      unique_file_ids,
      # TESTING PURPOSE ONLY
      # c("3cef5b6aa0f5772485bb710f71e69456___1.h5ad",
      #   "cd2caa6de850f73af4ca78a2ea307dd4___1.h5ad")
      get_unique_file_ids(cell_metadata) 
      # |> head(2)
      ,
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly", "HDF5Array"),
      deployment = "main"
    ),
    tar_target(
      file_id_cell_id_dict,
      create_file_id_cell_id_dict(cell_metadata, unique_file_ids),
      pattern = map(unique_file_ids),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly", "HDF5Array"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_5_minimal")
      )
    )
  )
  
}, script = paste0(store, "_target_script.R"), ask = FALSE)


job::job({
  
  tar_make(
    script = paste0(store, "_target_script.R"), 
    store = store, 
    reporter = "summary"
  )
  
})

file_id_cell_id_dict = tar_read(file_id_cell_id_dict, store = store)
file_id_cell_id_dict |> arrow::write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/file_id_cell_id_dict_v1_0_0.parquet",
                                             compression = "zstd")
rm(file_id_cell_id_dict)
gc()

