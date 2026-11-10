library(targets)
store = "/vast/scratch/users/shen.m/cellnexus_dataset_cell_map_Jul_2024_v1_2_1_target_store"
tar_script({
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(targets)
  library(tarchetypes)
  library(crew)
  library(crew.cluster)
  
  tar_option_set(
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue", 
    cue = tar_cue(mode = "never"),
    
    workspace_on_error = TRUE,
    controller = crew_controller_group(
      list(
        crew_controller_slurm(
          name = "elastic",
          workers = 300,
          tasks_max = 20,
          seconds_idle = 30,
          crashes_error = 10,
          options_cluster = crew_options_slurm(
            memory_gigabytes_required = c(25, 35, 40, 80, 160),
            cpus_per_task = c(2, 2, 5, 10, 20),
            time_minutes = c(30, 30, 30, 60*4, 60*24),
            verbose = T
          )
        )
      )
    ), 
    trust_object_timestamps = TRUE
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
    tar_target(cell_metadata , "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_5_0_mengyuan.parquet",
               deployment = "main"),
    tar_target(
      unique_file_ids,
      # TESTING PURPOSE ONLY
      # c("3cef5b6aa0f5772485bb710f71e69456___1.h5ad",
      #   "cd2caa6de850f73af4ca78a2ea307dd4___1.h5ad")
      get_unique_file_ids(cell_metadata) 
      # |> head(2)
      ,
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly", "HDF5Array")
      # resources = tar_resources(
      #   crew = tar_resources_crew(controller = "elastic")
      # )
    ),
    tar_target(
      file_id_cell_id_dict,
      create_file_id_cell_id_dict(cell_metadata, unique_file_ids),
      pattern = map(unique_file_ids),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly", "HDF5Array")
      # resources = tar_resources(
      #   crew = tar_resources_crew(controller = "elastic")
      # )
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
file_id_cell_id_dict |> arrow::write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/file_id_cell_id_dict_v1_1_1_Jul_2024.parquet",
                                     compression = "zstd")
rm(file_id_cell_id_dict)
gc()

