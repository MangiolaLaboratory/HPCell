#' @importFrom glue glue
#' @import targets
#'
#' @export
#' 
run_targets_pipeline <- function(
    input_data = file_path, 
    store =  store, 
    input_reference = input_reference_path,
    tissue,
    computing_resources = crew_controller_local(workers = 1), 
    debug_step = "annotation_label_transfer_tbl",
    filter_input = "TRUE", 
    RNA_assay_name = "originalexp"
    
){
  # Save inputs for passing to targets pipeline 
  input_data |> saveRDS("input_file.rds")
  input_reference |> saveRDS("input_reference.rds")
  tissue |> saveRDS("tissue.rds")
  computing_resources |> saveRDS("temp_computing_resources.rds")
  filtered |> saveRDS("filtered.rds")
  # Write pipeline to a file
  tar_script({
    
    library(targets)
    library(tarchetypes)
    library(crew)
    library(crew.cluster)
    
    computing_resources = readRDS("temp_computing_resources.rds")
    
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
      memory = "transient",
      garbage_collection = TRUE,
      #trust_object_timestamps = TRUE,
      storage = "worker", 
      retrieval = "worker", 
      #error = "continue",         
      format = "qs", 
      debug = debug_step, # Set the target you want to debug.
      # cue = tar_cue(mode = "never") # Force skip non-debugging outdated targets.
      controller = computing_resources
    )
    
    #-----------------------#
    # Future SLURM
    #-----------------------#
    
    # library(future)
    # library("future.batchtools")
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
    
    # small_slurm = 
    #   tar_resources(
    #     future = tar_resources_future(
    #       plan = tweak(
    #         batchtools_slurm,
    #         template = "dev/slurm_batchtools.tmpl",
    #         resources = list(
    #           ncpus = 2,
    #           memory = 40000,
    #           walltime = 172800
    #         )
    #       )
    #     )
    #   )
    # 
    # big_slurm = 
    #   tar_resources(
    #     future = tar_resources_future(
    #       plan = tweak(
    #         batchtools_slurm,
    #         template = "dev/slurm_batchtools.tmpl",
    #         resources = list(
    #           ncpus = 19,
    #           memory = 6000,
    #           walltime = 172800
    #         )
    #       )
    #     )
    #   )
    
    target_list = list(
      tar_target(file, "input_file.rds", format = "rds"), 
      tar_target(read_file, readRDS("input_file.rds")),
      tar_target(reference_file, "input_reference.rds", format = "rds"), 
      tar_target(read_reference_file, readRDS("input_reference.rds")), 
      tar_target(tissue_file, readRDS("tissue.rds")), 
      tar_target(filtered_file, readRDS("filtered.rds")))
    
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
      tar_target(filtered, filtered_file, deployment = "main"),
      tar_target(tissue, tissue_file, deployment = "main"),
      tar_target(reference_label_coarse, reference_label_coarse_id(tissue), deployment = "main"), 
      tar_target(reference_label_fine, reference_label_fine_id(tissue), deployment = "main"), 
      # Reading input files
      tar_target(input_read, readRDS(read_file),
                 pattern = map(read_file),
                 iteration = "list", deployment = "main"),
      tar_target(input_read_RNA_assay, add_RNA_assay(input_read, RNA_assay_name), 
                 pattern = map(input_read), 
                 iteration = "list"),
      tar_target(reference_read, readRDS(read_reference_file),
                 deployment = "main"),
      
      # Identifying empty droplets
      tar_target(empty_droplets_tbl,
                 empty_droplet_id(input_read_RNA_assay, filtered),
                 pattern = map(input_read_RNA_assay),
                 iteration = "list"),
      
      # Cell cycle scoring
      tar_target(cell_cycle_score_tbl, cell_cycle_scoring(input_read_RNA_assay,
                                                          empty_droplets_tbl),
                 pattern = map(input_read_RNA_assay,
                               empty_droplets_tbl),
                 iteration = "list"),
      
      # Annotation label transfer
      tar_target(annotation_label_transfer_tbl,
                 annotation_label_transfer(input_read_RNA_assay,
                                           reference_read,
                                           empty_droplets_tbl),
                 pattern = map(input_read_RNA_assay,
                               empty_droplets_tbl),
                 iteration = "list"),
      
      # Alive identification
      tar_target(alive_identification_tbl, alive_identification(input_read_RNA_assay,
                                                                empty_droplets_tbl,
                                                                annotation_label_transfer_tbl),
                 pattern = map(input_read_RNA_assay,
                               empty_droplets_tbl,
                               annotation_label_transfer_tbl),
                 iteration = "list"),
      
      # Doublet identification
      tar_target(doublet_identification_tbl, doublet_identification(input_read_RNA_assay,
                                                                    empty_droplets_tbl,
                                                                    alive_identification_tbl,
                                                                    annotation_label_transfer_tbl,
                                                                    reference_label_fine),
                 pattern = map(input_read_RNA_assay,
                               empty_droplets_tbl,
                               alive_identification_tbl,
                               annotation_label_transfer_tbl),
                 iteration = "list"),
      
      # Non-batch variation removal
      tar_target(non_batch_variation_removal_S, non_batch_variation_removal(input_read_RNA_assay,
                                                                            empty_droplets_tbl,
                                                                            alive_identification_tbl,
                                                                            cell_cycle_score_tbl),
                 pattern = map(input_read_RNA_assay,
                               empty_droplets_tbl,
                               alive_identification_tbl,
                               cell_cycle_score_tbl),
                 iteration = "list"),
      
      # Pre-processing output
      tar_target(preprocessing_output_S, preprocessing_output(tissue,
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
                 iteration = "list"),
      
      # pseudobulk preprocessing
      tar_target(pseudobulk_preprocessing_SE, pseudobulk_preprocessing(reference_label_fine,
                                                                       preprocessing_output_S))
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
    store = store, 
    callr_function = NULL
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
