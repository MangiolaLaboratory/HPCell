
#' @importFrom glue glue
#' @import targets
#'
#' @export
#' 
run_targets_pipeline <- function(input_data, store =  tempfile(tmpdir = "."), input_reference, tissue, reference_label){
 
  # Save inputs for passing to targets pipeline 
  input_data |> saveRDS("input_file.rds")
  input_reference |> saveRDS("input_reference.rds")
  tissue |> save(tissue, file = "tissue.RData")
  reference_label |> save(reference_label, file = "reference_label.RData")
  # Write pipeline to a file
  tar_script({
    library(targets)
    library(tarchetypes)
    # tar_source(files = "R")
    #library(Seurat)
    #library(tidyseurat)
    #source("/home/users/allstaff/si.j/jascap/dev/targets_jascap/R/Packages.R")
    # source("/home/users/allstaff/si.j/jascap/dev/targets_jascap/R/Functions.R")
    #-----------------------#
    # Packages
    #-----------------------#
    tar_option_set(
      packages = c(
        "jascap",
        "readr",
        "dplyr",
        "tidyr",
        "ggplot2",
        "purrr",
        "Seurat",
        "tidyseurat",
        "glue",
        "scater",
        "DropletUtils",
        "EnsDb.Hsapiens.v86",
        "here",
        "stringr",
        "readr",
        "rlang",
        "scuttle",
        "scDblFinder",
        "ggupset",
        "tidySummarizedExperiment",
        "tidyverse",
        "broom",
        "tarchetypes",
        "SeuratObject",
        "SingleCellExperiment", 
        "SingleR", 
        "celldex", 
        "tidySingleCellExperiment", 
        "tibble", 
        "magrittr",
        "qs", 
        "S4Vectors"
      ),
      storage = "worker", 
      retrieval = "worker", 
      #error = "continue",         
      format = "qs"
    )
    
    #-----------------------#
    # Future SLURM
    #-----------------------#
    
    library(future)
    library("future.batchtools")
    # slurm <- 
    #     `batchtools_slurm` |>
    #     future::tweak( template = glue("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/third_party_sofware/slurm_batchtools.tmpl"),
    #                                  resources=list(
    #                                      ncpus = 20,
    #                                      memory = 6000,
    #                                      walltime = 172800
    #                                  )
    #     )
    # plan(slurm)
    
    small_slurm = 
      tar_resources(
        future = tar_resources_future(
          plan = tweak(
            batchtools_slurm,
            template = "dev/slurm_batchtools.tmpl",
            resources = list(
              ncpus = 2,
              memory = 40000,
              walltime = 172800
            )
          )
        )
      )
    
    big_slurm = 
      tar_resources(
        future = tar_resources_future(
          plan = tweak(
            batchtools_slurm,
            template = "dev/slurm_batchtools.tmpl",
            resources = list(
              ncpus = 19,
              memory = 6000,
              walltime = 172800
            )
          )
        )
      )
    
    target_list = list(
      tar_target(file, "input_file.rds", format = "rds"), 
      tar_target(read_file, readRDS("input_file.rds")),
      tar_target(reference_file, "input_reference.rds", format = "rds"), 
      tar_target(read_reference_file, readRDS("input_reference.rds")), 
      tar_target(tissue_load, load("tissue.RData")),
      tar_target(reference_label_load, load("reference_label.RData")))
    
    #-----------------------#
    # Pipeline
    #-----------------------#
    target_list|> c(list(
      
      # Define input files
      # tarchetypes::tar_files(name= input_track, 
      #                        read_file, 
      #                        deployment = "main"),
      # tarchetypes::tar_files(name= reference_track,
      #                        read_reference_file, 
      #                        deployment = "main"),
      tar_target(filtered, "filtered", deployment = "main"),
      tar_target(tissue, tissue_load, deployment = "main"),
      tar_target(reference_label, reference_label_load, deployment = "main"),
      tar_target(reference_label_identification, reference_label_id(tissue, reference_label), deployment = "main"), 
      # Reading input files
      tar_target(input_read, readRDS(read_file),
                 pattern = map(read_file),
                 iteration = "list", deployment = "main"),
      tar_target(reference_read, readRDS(read_reference_file),
                 deployment = "main"),
      
      # Identifying empty droplets
      tar_target(empty_droplets_tbl,
                 empty_droplet_id(input_read, filtered),
                 pattern = map(input_read),
                 iteration = "list", resources = small_slurm),
      
      # Cell cycle scoring
      tar_target(cell_cycle_score_tbl, cell_cycle_scoring(input_read,
                                                      empty_droplets_tbl),
                 pattern = map(input_read,
                               empty_droplets_tbl),
                 iteration = "list", resources = small_slurm),
      
      # Annotation label transfer
      tar_target(annotation_label_transfer_tbl,
                 annotation_label_transfer(input_read,
                                           reference_read,
                                           empty_droplets_tbl),
                 pattern = map(input_read,
                               empty_droplets_tbl),
                 iteration = "list", resources = small_slurm),
      
      # Alive identification
      tar_target(alive_identification_tbl, alive_identification(input_read,
                                                            empty_droplets_tbl,
                                                            annotation_label_transfer_tbl),
                 pattern = map(input_read,
                               empty_droplets_tbl,
                               annotation_label_transfer_tbl),
                 iteration = "list", resources = small_slurm),
      
      # Doublet identification
      tar_target(doublet_identification_tbl, doublet_identification(input_read,
                                                                empty_droplets_tbl,
                                                                alive_identification_tbl,
                                                                annotation_label_transfer_tbl,
                                                                reference_label_identification),
                 pattern = map(input_read,
                               empty_droplets_tbl,
                               alive_identification_tbl,
                               annotation_label_transfer_tbl),
                 iteration = "list", resources = small_slurm),
      
      # Non-batch variation removal
      tar_target(non_batch_variation_removal_S, non_batch_variation_removal(input_read,
                                                                          empty_droplets_tbl,
                                                                          alive_identification_tbl,
                                                                          cell_cycle_score_tbl),
                 pattern = map(input_read,
                               empty_droplets_tbl,
                               alive_identification_tbl,
                               cell_cycle_score_tbl),
                 iteration = "list", resources = small_slurm),
      
      # Pre-processing output
      tar_target(preprocessing_output_S, preprocessing_output("none",
                                                            non_batch_variation_removal_S,
                                                            alive_identification_tbl,
                                                            cell_cycle_score_tbl,
                                                            annotation_label_transfer_tbl,
                                                            doublet_identification_tbl),
                 pattern = map(non_batch_variation_removal_S,
                               alive_identification_tbl,
                               cell_cycle_score_tbl,
                               annotation_label_transfer_tbl,
                               doublet_identification_tbl),
                 iteration = "list", resources = small_slurm),
      
      # pseudobulk preprocessing
      tar_target(pseudobulk_preprocessing_SE, pseudobulk_preprocessing(reference_label_identification,
                                                                    preprocessing_output_S),
                 resources = small_slurm)
    ))
    
  }, script = glue("{store}.R"), ask = FALSE)
  
  #Running targets 
  # input_files<- c("CB150T04X__batch14.rds","CB291T01X__batch8.rds")
  # run_targets <- function(input_files){
  #   tar_make(
  #     script = glue("{store}.R"),
  #     store = store
  #   )
  # }
  # run_targets(input_files)
  tar_make(
    script = glue("{store}.R"),
    store = store
  )
  # tar_make_future(
  #   script = glue("{store}.R"),
  #   store = store, 
  #   workers = 200, 
  #   garbage_collection = TRUE
  # )
  
  message(glue("HPCell says: you can read your output executing tar_read(preprocessing_output_S, store = \"{store}\") "))
  
  tar_read(preprocessing_output_S, store = store)
  
}

## my_results = run_targets_pipeline(..)
