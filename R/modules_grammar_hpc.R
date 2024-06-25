#' Run Targets Pipeline for HPCell
#'
#' @description
#' This function sets up and executes a `targets` pipeline for HPCell. It saves input data and configurations,
#' writes a pipeline script, and runs the pipeline using the 'targets' package.
#'
#' @param input_data Input data for the pipeline.
#' @param store Directory path for storing the pipeline files.
#' @param input_reference Optional reference data.
#' @param tissue Tissue type for the analysis.
#' @param computing_resources Configuration for computing resources.
#' @param debug_step Optional step for debugging.
#' @param filter_empty_droplets Flag to indicate if input filtering is needed.
#' @param RNA_assay_name Name of the RNA assay.
#' @param sample_column Column name for sample identification.
#' @param cell_type_annotation_column Column name for cell type annotation in input data
#'
#' @return The output of the `targets` pipeline, typically a pre-processed data set.
#'
#' @importFrom glue glue
#' @importFrom targets tar_script
#' @importFrom magrittr set_names
#' @import crew.cluster
#' @import tarchetypes
#' @import targets
#' @import broom
#' @import ggplot2
#' @import ggupset
#' @import here
#' @import qs
#' @import crew
#' @importFrom future tweak
#' @import crew
#' @import crew.cluster
#' @export
initialise_hpc <- function(input_data,
                           store =  targets::tar_config_get("store"),
                           computing_resources = crew_controller_local(workers = 1),
                           debug_step = NULL,
                           RNA_assay_name = "RNA") {
  # Capture all arguments including defaults
  args_list <- as.list(environment())
  
  # if simple names are not set, use integers
  if(input_data |> names() |> is.null())
    input_data |> set_names(seq_len(length(input_data)))
  
  input_data |> names() |> saveRDS("sample_names.rds")
  

  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  # Add pipeline step
  args_list$target_chunk = 
    substitute({
      
      target_list = 
        target_list |> c(list(
          tar_target(sample_name, readRDS("sample_names.rds"))
        ))
      
    })
  
  
  list(initialisation = args_list) |>
    add_class("HPCell")
}


target_chunk_remove_empty_DropletUtils = 
  substitute({
    target_list = 
      target_list |> c(list(
        tar_target(my_total_RNA_count_check, readRDS("total_RNA_count_check.rds")), 
        tar_target(
          empty_droplets_tbl,
          empty_droplet_id(input_read, my_total_RNA_count_check),
          pattern = map(input_read),
          iteration = "list"
        )
      ))
    
  })

target_chunk_undefined_remove_empty_DropletUtils = 
  substitute({
    target_list = 
      target_list |> c(list(
        tar_target(
          empty_droplets_tbl, 
          input_read |> as_tibble() |> select(.cell) |> mutate(empty_droplet = FALSE), 
          pattern = map(input_read), 
          iteration = "list",
          packages = c("dplyr","tidySingleCellExperiment", "tidyseurat")
          
        )
      ))
    
  })


# Define the generic function
#' @export
remove_empty_DropletUtils <- function(input_data, total_RNA_count_check = NULL, ...) {
  UseMethod("remove_empty_DropletUtils")
}

#' @export
remove_empty_DropletUtils.Seurat = function(input_data, total_RNA_count_check = NULL, ...) {
  # Capture all arguments including defaults
  args_list <- as.list(environment())
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  list(initialisation = list(input_data = input_data)) |>
    add_class("HPCell") |>
    remove_empty_DropletUtils()
  
}


#' @export
remove_empty_DropletUtils.HPCell = function(input_hpc, total_RNA_count_check = NULL, ...) {
  
  # Capture all arguments including defaults
  args_list <- as.list(environment())[-1]
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  # Save parameter
  total_RNA_count_check |> saveRDS("total_RNA_count_check.rds")  
  
  # Add pipeline step
  input_hpc |>
    c(list(remove_empty_DropletUtils = args_list)) |>
    add_class("HPCell")
  
  
}


target_chunk_undefined_remove_dead_scuttle = 
  substitute({
    target_list = 
      target_list |> c(list(
        tar_target(
          alive_identification_tbl, 
          input_read |> as_tibble() |> dplyr::select(.cell) |> mutate(alive = TRUE), 
          pattern = map(input_read), 
          iteration = "list",
          packages = c("dplyr","tidySingleCellExperiment", "tidyseurat")
          
        )
      ))
    
  })




# Define the generic function
#' @export
remove_dead_scuttle <- function(input_data, group_by = NULL) {
  UseMethod("remove_dead_scuttle")
}

#' @export
remove_dead_scuttle.HPCell = function(input_hpc, group_by = NULL) {
  
  # Capture all arguments including defaults
  args_list <- as.list(environment())[-1]
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  group_by |> saveRDS("temp_group_by.rds")
  
  args_list$target_chunk = 
    substitute({
      target_list = 
        target_list |> c(list(
          tar_target(grouping_column, readRDS("temp_group_by.rds")),
          tar_target(
            alive_identification_tbl, 
            alive_identification(input_read,
                                 empty_droplets_tbl,
                                 annotation_label_transfer_tbl,
                                 grouping_column),
            pattern = map(input_read,
                          empty_droplets_tbl,
                          annotation_label_transfer_tbl),
            iteration = "list"
          )
        ))
      
    })
  
  input_hpc |>
    c(list(remove_dead_scuttle = args_list))  |>
    add_class("HPCell")
  
}

target_chunk_undefined_score_cell_cycle_seurat = 
  substitute({
    target_list = 
      target_list |> c(list(
        tar_target(cell_cycle_score_tbl,
                   NULL,
                   pattern = map(input_read),
                   iteration = "list", 
                   description = "setup"
        )
      ))
    
  })

# Define the generic function
#' @export
score_cell_cycle_seurat <- function(input_data, ...) {
  UseMethod("score_cell_cycle_seurat")
}

#' @export
score_cell_cycle_seurat.HPCell = function(input_hpc) {
  # Capture all arguments including defaults
  args_list <- as.list(environment())[-1]
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  # Add pipeline step
  args_list$target_chunk = 
    substitute({
      
      target_list = 
        target_list |> c(list(
          tar_target(
            cell_cycle_score_tbl,
            cell_cycle_scoring(input_read,
                               empty_droplets_tbl),
            pattern = map(input_read,
                          empty_droplets_tbl),
            iteration = "list"
          )
        ))
      
    })
  
  input_hpc |>
    c(list(score_cell_cycle_seurat = args_list)) |>
    add_class("HPCell")
  
}




target_chunk_undefined_remove_doublets_scDblFinder = 
  substitute({
    target_list = 
      target_list |> c(list(
        tar_target(
          doublet_identification_tbl, 
          input_read |> as_tibble() |> dplyr::select(.cell) |> mutate(scDblFinder.class=="singlet"), 
          pattern = map(input_read), 
          iteration = "list",
          packages = c("dplyr","tidySingleCellExperiment", "tidyseurat")
          
        )
      ))
    
  })




# Define the generic function
#' @export
remove_doublets_scDblFinder <- function(input_hpc) {
  UseMethod("remove_doublets_scDblFinder")
}

#' @export
remove_doublets_scDblFinder.HPCell = function(input_hpc) {
  # Capture all arguments including defaults
  args_list <- as.list(environment())[-1]
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  args_list$target_chunk = 
    substitute({
      target_list = 
        target_list |> c(list(
          tar_target(doublet_identification_tbl, doublet_identification(input_read,
                                                                        empty_droplets_tbl,
                                                                        alive_identification_tbl),
                     pattern = map(input_read,
                                   empty_droplets_tbl,
                                   alive_identification_tbl
                     ),
                     iteration = "list"
          )
        ))
      
    })
  
  input_hpc |>
    c(list(remove_doublets_scDblFinder = args_list)) |>
    add_class("HPCell")
  
}

target_chunk_annotate_cell_type = 
  substitute({
    target_list = 
      target_list |> c(list(
        tar_target(reference_read, readRDS("input_reference.rds")), 
        tar_target(annotation_label_transfer_tbl,
                   annotation_label_transfer(input_read,
                                             empty_droplets_tbl,
                                             reference_read),
                   pattern = map(input_read, empty_droplets_tbl),
                   iteration = "list"
        )
        
      ))
    
  })

target_chunk_undefined_annotate_cell_type = 
  substitute({
    target_list = 
      target_list |> c(list(
        tar_target(annotation_label_transfer_tbl,
                   NULL,
                   pattern = map(input_read),
                   iteration = "list"
        )
      ))
    
  })


# Define the generic function
#' @export
annotate_cell_type <- function(input_data, ...) {
  UseMethod("annotate_cell_type")
}

#' @export
annotate_cell_type.HPCell = function(input_hpc, azimuth_reference = NULL) {
  # Capture all arguments including defaults
  args_list <- as.list(environment())[-1]
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  input_hpc |>
    c(list(annotate_cell_type = args_list)) |>
    add_class("HPCell")
  
}


target_chunk_undefined_normalise_abundance_seurat_SCT = 
  substitute({
    target_list = 
      target_list |> c(list(
        tar_target( 
          non_batch_variation_removal_S, 
          NULL, 
          pattern = map(input_read), 
          iteration = "list"
          
        )
      ))
    
  })


# Define the generic function
#' @export
normalise_abundance_seurat_SCT <- function(input_data, ...) {
  UseMethod("normalise_abundance_seurat_SCT")
}

#' @export
normalise_abundance_seurat_SCT.HPCell = function(input_hpc, factors_to_regress = NULL) {
  # Capture all arguments including defaults
  args_list <- as.list(environment())[-1]
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  factors_to_regress |> saveRDS("factors_to_regress.rds")
  
  args_list$target_chunk = 
    substitute({
      target_list = 
        target_list |> c(list(
          tar_target(my_factors_to_regress, readRDS("factors_to_regress.rds")), 
          tar_target(non_batch_variation_removal_S, non_batch_variation_removal(input_read,
                                                                                empty_droplets_tbl,
                                                                                alive_identification_tbl,
                                                                                cell_cycle_score_tbl,
                                                                                factors_to_regress = my_factors_to_regress),
                     pattern = map(input_read,
                                   empty_droplets_tbl,
                                   alive_identification_tbl,
                                   cell_cycle_score_tbl),
                     iteration = "list"
          )
        ))
      
    })
  
  input_hpc |>
    c(list(normalise_abundance_seurat_SCT = args_list)) |>
    add_class("HPCell")
  
}


# Define the generic function
#' @export
calculate_pseudobulk <- function(input_data, ...) {
  UseMethod("calculate_pseudobulk")
}

#' @export
calculate_pseudobulk.HPCell = function(input_hpc, group_by = NULL) {
  # Capture all arguments including defaults
  args_list <- as.list(environment())[-1]
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  group_by |> saveRDS("pseudobulk_group_by.rds")
  
  args_list$target_chunk = 
    substitute({
      target_list = 
        target_list |> c(list(
          tar_target(pseudobulk_group_by, readRDS("pseudobulk_group_by.rds")), 
          tar_target(create_pseudobulk_sample, create_pseudobulk(
            preprocessing_output_S,  
            sample_names,
            x = pseudobulk_group_by
          ), 
                     pattern = map(preprocessing_output_S), 
                     iteration = "list"),
          tar_target(
            pseudobulk_merge_all_samples, 
            create_pseudobulk_sample |> pseudobulk_merge()
          )
        ))
      
    })
  
  input_hpc |>
    c(list(calculate_pseudobulk = args_list)) |>
    add_class("HPCell")
  
}




# pseudobulk preprocessing for each sample 





# Define the generic function
#' @export
evaluate_hpc <- function(input_data) {
  UseMethod("evaluate_hpc")
}

#' @importFrom glue glue
#' @export
evaluate_hpc.HPCell = function(input_hpc) {
# 
#   # Fix GCHECKS
#   read_file <- NULL
#   reference_file <- NULL
#   tissue_file <- NULL
#   filtered_file <- NULL
#   sample_column_file <- NULL
#   cell_type_annotation_column_file <- NULL
#   reference_label_coarse <- NULL
#   reference_label_fine <- NULL
#   input_read <- NULL
#   unique_tissues <- NULL
#   reference_read <- NULL
#   empty_droplets_tbl <- NULL
#   cell_cycle_score_tbl <- NULL
#   annotation_label_transfer_tbl <- NULL
#   alive_identification_tbl <- NULL
#   doublet_identification_tbl <- NULL
#   non_batch_variation_removal_S <- NULL
#   preprocessing_output_S <- NULL
#   create_pseudobulk_sample <- NULL
#   sampleName <- NULL
#   cellAnno <- NULL
#   pseudobulk_merge_all_samples <- NULL
#   calc_UMAP_dbl_report <- NULL
#   variable_gene_list <- NULL
#   tar_render <- NULL
#   empty_droplets_report <- NULL
#   doublet_identification_report <- NULL
#   Technical_variation_report <- NULL
#   pseudobulk_processing_report <- NULL
  
  #sample_column = enquo(input_hpc$initialisation$sample_column)
  # cell_type_annotation_column = enquo(cell_type_annotation_column)
  
  # Save inputs for passing to targets pipeline
  # input_data |> CHANGE_ASSAY |> saveRDS("input_file.rds")
  
  dir.create(input_hpc$initialisation$store, showWarnings = FALSE, recursive = TRUE)
  data_file_names = glue("{input_hpc$initialisation$store}/{names(input_hpc$initialisation$input_data)}.rds")
    map2(
      input_hpc$initialisation$input_data,
      data_file_names,
      ~ .x |> saveRDS(.y)
    )
  
    list(data_file_names) |> saveRDS("input_file.rds")
  
  input_hpc$initialisation$computing_resources |> saveRDS("temp_computing_resources.rds")
  #sample_column |> saveRDS("sample_column.rds")
  input_hpc$initialisation$debug_step |> saveRDS("temp_debug_step.rds")
  
  # Write pipeline to a file
  tar_script({
    library(targets)
    library(tarchetypes)
    library(crew)
    library(crew.cluster)
    
    tar_option_set(
      memory = "transient",
      garbage_collection = TRUE,
      storage = "worker",
      retrieval = "worker",
      #error = "continue",
      format = "qs",
      debug = readRDS("temp_debug_step.rds"), # Set the target you want to debug.
      # cue = tar_cue(mode = "never") # Force skip non-debugging outdated targets.
      controller = readRDS("temp_computing_resources.rds"), 
      packages = c("HPCell")
    )
    
    target_list = list(
      tar_target(file, "input_file.rds", format = "rds"),
      tar_target(read_file, readRDS(file), iteration = "list")
    )
    
    target_list = 
      target_list |> c(list(
        
        # Reading input files
        tar_target(
          input_read,
          readRDS(read_file),
          pattern = map(read_file),
          iteration = "list"
        )
      ))
    
  }, script = glue("{input_hpc$initialisation$store}.R"), ask = FALSE)
  
  # Set sample names
  tar_script_append2(
    input_hpc$initialisation$target_chunk,
    script = glue("{input_hpc$initialisation$store}.R")
  )
  
  #-----------------------#
  # Empty droplets
  #-----------------------#
  
  if("remove_empty_DropletUtils" %in% names(input_hpc))
    tar_script_append2(
      target_chunk_remove_empty_DropletUtils,
      script = glue("{input_hpc$initialisation$store}.R")
    )
  else
    tar_script_append2(
      target_chunk_undefined_remove_empty_DropletUtils,
      script = glue("{input_hpc$initialisation$store}.R")
    )
  
  #-----------------------#
  # Annotate cell type
  #-----------------------#
  
  if(
    "annotate_cell_types" %in% names(input_hpc) |
    ( "remove_dead_scuttle" %in% names(input_hpc) & !is.null(input_hpc$remove_dead_scuttle$group_by))
  )
    tar_script_append2(
      target_chunk_annotate_cell_type,
      script = glue("{input_hpc$initialisation$store}.R")
    )
  else
    tar_script_append2(
      target_chunk_undefined_annotate_cell_type,
      script = glue("{input_hpc$initialisation$store}.R")
    )
  
  #-----------------------#
  # Remove dead
  #-----------------------#
  
  if("remove_dead_scuttle" %in% names(input_hpc))
    tar_script_append2(
      input_hpc$remove_dead_scuttle$target_chunk,
      script = glue("{input_hpc$initialisation$store}.R")
    )
  else
    tar_script_append2(
      target_chunk_undefined_remove_dead_scuttle,
      script = glue("{input_hpc$initialisation$store}.R")
    )
  
  
  #-----------------------#
  # score cell cycle
  #-----------------------#
  if("score_cell_cycle_seurat" %in% names(input_hpc))
    tar_script_append2(
      input_hpc$score_cell_cycle_seurat$target_chunk,
      script = glue("{input_hpc$initialisation$store}.R")
    )
  else
    tar_script_append2(
      target_chunk_undefined_score_cell_cycle_seurat,
      script = glue("{input_hpc$initialisation$store}.R")
    )
  
  #-----------------------#
  # Doublets
  #-----------------------#
  
  if("remove_doublets_scDblFinder" %in% names(input_hpc))
    tar_script_append2(
      input_hpc$remove_doublets_scDblFinder$target_chunk,
      script = glue("{input_hpc$initialisation$store}.R")
    )
  else
    tar_script_append2(
      target_chunk_undefined_remove_doublets_scDblFinder,
      script = glue("{input_hpc$initialisation$store}.R")
    )
  
  #-----------------------#
  # SCT
  #-----------------------#
  
  if("normalise_abundance_seurat_SCT" %in% names(input_hpc))
    tar_script_append2(
      input_hpc$normalise_abundance_seurat_SCT$target_chunk,
      script = glue("{input_hpc$initialisation$store}.R")
    )
  else
    tar_script_append2(
      target_chunk_undefined_normalise_abundance_seurat_SCT,
      script = glue("{input_hpc$initialisation$store}.R")
    )
  
  
  
  #-----------------------#
  # Create single cell output
  #-----------------------#
  
  # Pre-processing output
  
  
  # Pre-processing output
  tar_script_append({
    target_list = 
      target_list |> c(list(
        
        tar_target(preprocessing_output_S, preprocessing_output(input_read,
                                                                empty_droplets_tbl,
                                                                non_batch_variation_removal_S,
                                                                alive_identification_tbl,
                                                                cell_cycle_score_tbl,
                                                                annotation_label_transfer_tbl,
                                                                doublet_identification_tbl),
                   pattern = map(input_read,
                                 empty_droplets_tbl,
                                 non_batch_variation_removal_S,
                                 alive_identification_tbl,
                                 cell_cycle_score_tbl,
                                 annotation_label_transfer_tbl,
                                 doublet_identification_tbl),
                   iteration = "list")
      ))
  }, script = glue("{input_hpc$initialisation$store}.R"))
  
  
  #-----------------------#
  # Create pseudobulk 
  #-----------------------#
  
  if("calculate_pseudobulk" %in% names(input_hpc))
    tar_script_append2(
      input_hpc$calculate_pseudobulk$target_chunk,
      script = glue("{input_hpc$initialisation$store}.R")
    )

  
  #-----------------------#
  # Close pipeline
  #-----------------------#
  
  # Call final list
  tar_script_append({
    target_list
  }, script = glue("{input_hpc$initialisation$store}.R"))
  
  if(readRDS("temp_debug_step.rds") |> is.null())
    my_callr_function =  callr::r
  else
    my_callr_function =  NULL
  
  tar_make(
    callr_function = my_callr_function,
    reporter = "verbose_positives",
    script = glue("{input_hpc$initialisation$store}.R"),
    store = input_hpc$initialisation$store
  )
  
  # Example usage:
  c(
    "input_file.rds",
    "temp_computing_resources.rds",
    "temp_debug_step.rds",
    "sample_names.rds",
    "total_RNA_count_check.rds",
    "temp_group_by.rds",
    "factors_to_regress.rds",
    "pseudobulk_group_by.rds"
  ) |> 
    remove_files_safely()
  
  return(
    tar_meta(preprocessing_output_S, store = glue("{input_hpc$initialisation$store}"))
  )
}

#' @importFrom methods show
#' @export
print.HPCell <- function(x, ...){
  
  x |>
    evaluate_hpc() |> 
    print()
}
