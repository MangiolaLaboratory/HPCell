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
                           tier = NULL,
                           debug_step = NULL,
                           RNA_assay_name = "RNA") {
  # Capture all arguments including defaults
  args_list <- as.list(environment())
  
  # if simple names are not set, use integers
  if(input_data |> names() |> is.null())
    input_data |> set_names(seq_len(length(input_data)))
  
  input_data |> names() |> saveRDS("sample_names.rds")
  #cell_count |> saveRDS("cell_count.rds")

  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  # Add pipeline step
  args_list$target_chunk = 
    substitute({
      
      target_list = 
        target_list |> c(list(
          tar_target(sample_names, readRDS("sample_names.rds"))
        ))
      
    })
  
  
  list(initialisation = args_list) |>
    add_class("HPCell")
}





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
  args_list$target_chunk = function(){
    
    append_chunk_fix(
      { tar_target(my_total_RNA_count_check, readRDS("total_RNA_count_check.rds")) }, 
      script = glue("{input_hpc$initialisation$store}.R")
    )
    
    append_chunk_tiers(
      { tar_target(
        empty_droplets_tbl_TIER_PLACEHOLDER,
        empty_droplet_id(input_read, my_total_RNA_count_check),
        pattern = slice(input_read, index  = SLICE_PLACEHOLDER ),
        iteration = "list", 
        resources = RESOURCE_PLACEHOLDER
      ) }, 
      tiers = input_hpc$initialisation$tier,
      script = glue("{input_hpc$initialisation$store}.R")
    )
  
  }
  
  # Add pipeline step
  input_hpc |>
    c(list(remove_empty_DropletUtils = args_list)) |>
    add_class("HPCell")
  
  
}

target_chunk_undefined_remove_empty_DropletUtils = function(){
  append_chunk_tiers(
    { tar_target(
      empty_droplets_tbl_TIER_PLACEHOLDER,
      input_read |> as_tibble() |> select(.cell) |> mutate(empty_droplet = FALSE),
      pattern = slice(input_read, index  = SLICE_PLACEHOLDER ),
      iteration = "list", 
      resources = RESOURCE_PLACEHOLDER
    ) }, 
    tiers = input_hpc$initialisation$tier,
    script = glue("{input_hpc$initialisation$store}.R")
  )
}



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
  
  args_list$target_chunk = function(){
    
    append_chunk_fix(
      { tar_target(grouping_column, readRDS("temp_group_by.rds")) }, 
      script = glue("{input_hpc$initialisation$store}.R")
    )
    
    append_chunk_tiers(
      { tar_target(
        alive_identification_tbl_TIER_PLACEHOLDER, 
        alive_identification(input_read,
                             empty_droplets_tbl_TIER_PLACEHOLDER,
                             annotation_label_transfer_tbl_TIER_PLACEHOLDER,
                             grouping_column),
        pattern = map(slice(input_read, index  = SLICE_PLACEHOLDER ),
                      empty_droplets_tbl_TIER_PLACEHOLDER,
                      annotation_label_transfer_tbl_TIER_PLACEHOLDER),
        iteration = "list",
        resources = RESOURCE_PLACEHOLDER
      ) }, 
      tiers = input_hpc$initialisation$tier,
      script = glue("{input_hpc$initialisation$store}.R")
    )
    
  }
  
  input_hpc |>
    c(list(remove_dead_scuttle = args_list))  |>
    add_class("HPCell")
  
}


target_chunk_undefined_remove_dead_scuttle = function(input_hpc){
  append_chunk_tiers(
    { tar_target(
      alive_identification_tbl_TIER_PLACEHOLDER, 
      input_read |> as_tibble() |> dplyr::select(.cell) |> mutate(alive = TRUE), 
      pattern = slice(input_read, index  = SLICE_PLACEHOLDER ), 
      iteration = "list", 
      resources = RESOURCE_PLACEHOLDER
    ) }, 
    tiers = input_hpc$initialisation$tier,
    script = glue("{input_hpc$initialisation$store}.R")
  )
}



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
  args_list$target_chunk = function(){
    
    append_chunk_tiers(
      { tar_target(
        cell_cycle_score_tbl_TIER_PLACEHOLDER,
        cell_cycle_scoring(input_read,
                           empty_droplets_tbl_TIER_PLACEHOLDER),
        pattern = map(slice(input_read, index  = SLICE_PLACEHOLDER ),
                      empty_droplets_tbl_TIER_PLACEHOLDER),
        iteration = "list",
        resources = RESOURCE_PLACEHOLDER
      
      ) }, 
      tiers = input_hpc$initialisation$tier,
      script = glue("{input_hpc$initialisation$store}.R")
    )
    
  }
  
  input_hpc |>
    c(list(score_cell_cycle_seurat = args_list)) |>
    add_class("HPCell")
  
}

target_chunk_undefined_score_cell_cycle_seurat = function(input_hpc){
  append_chunk_tiers(
    { tar_target(
      cell_cycle_score_tbl_TIER_PLACEHOLDER,
      NULL,
      pattern = slice(input_read, index  = SLICE_PLACEHOLDER ),
      iteration = "list",
      resources = RESOURCE_PLACEHOLDER
    ) }, 
    tiers = input_hpc$initialisation$tier,
    script = glue("{input_hpc$initialisation$store}.R")
  )
}




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
  
  args_list$target_chunk = function(){
    
    append_chunk_tiers(
      { tar_target(
        doublet_identification_tbl_TIER_PLACEHOLDER, 
        doublet_identification(input_read, empty_droplets_tbl_TIER_PLACEHOLDER, alive_identification_tbl_TIER_PLACEHOLDER),
        pattern = map(slice(input_read, index  = SLICE_PLACEHOLDER ),
                      empty_droplets_tbl_TIER_PLACEHOLDER,
                      alive_identification_tbl_TIER_PLACEHOLDER
        ),
        iteration = "list",
        resources = RESOURCE_PLACEHOLDER
        
      ) }, 
      tiers = input_hpc$initialisation$tier,
      script = glue("{input_hpc$initialisation$store}.R")
    )
    
  }
  
  
  input_hpc |>
    c(list(remove_doublets_scDblFinder = args_list)) |>
    add_class("HPCell")
  
}

target_chunk_undefined_remove_doublets_scDblFinder = function(input_hpc){
  append_chunk_tiers(
    { tar_target(
      doublet_identification_tbl_TIER_PLACEHOLDER, 
      input_read |> as_tibble() |> dplyr::select(.cell) |> mutate(scDblFinder.class="singlet"), 
      pattern = slice(input_read, index  = SLICE_PLACEHOLDER ), 
      iteration = "list", 
      resources = RESOURCE_PLACEHOLDER
    ) }, 
    tiers = input_hpc$initialisation$tier,
    script = glue("{input_hpc$initialisation$store}.R")
  )
}

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
  
  args_list$target_chunk = function(){
    
    append_chunk_fix(
      {  tar_target(reference_read, readRDS("input_reference.rds")) }, 
      script = glue("{input_hpc$initialisation$store}.R")
    )
    
    append_chunk_tiers(
      { tar_target(annotation_label_transfer_tbl_TIER_PLACEHOLDER,
                   annotation_label_transfer(input_read,
                                             empty_droplets_tbl_TIER_PLACEHOLDER,
                                             reference_read),
                   pattern = map(slice(input_read, index  = SLICE_PLACEHOLDER ), empty_droplets_tbl_TIER_PLACEHOLDER),
                   iteration = "list",
                   resources = RESOURCE_PLACEHOLDER
        )}, 
      tiers = input_hpc$initialisation$tier,
      script = glue("{input_hpc$initialisation$store}.R")
    )
    
  }
  
  input_hpc |>
    c(list(annotate_cell_type = args_list)) |>
    add_class("HPCell")
  
}

target_chunk_undefined_annotate_cell_type = function(input_hpc){
  append_chunk_tiers(
    { tar_target(
      annotation_label_transfer_tbl_TIER_PLACEHOLDER, 
      NULL, 
      pattern = slice(input_read, index  = SLICE_PLACEHOLDER ), 
      iteration = "list", 
      resources = RESOURCE_PLACEHOLDER
    ) }, 
    tiers = input_hpc$initialisation$tier,
    script = glue("{input_hpc$initialisation$store}.R")
  )
}


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
  
  args_list$target_chunk = function(){
    
    append_chunk_fix(
      { tar_target(my_factors_to_regress, readRDS("factors_to_regress.rds")) }, 
      script = glue("{input_hpc$initialisation$store}.R")
    )
    
    append_chunk_tiers(
      { tar_target(
        non_batch_variation_removal_S_TIER_PLACEHOLDER, 
        non_batch_variation_removal(
              input_read,
             empty_droplets_tbl_TIER_PLACEHOLDER,
             alive_identification_tbl_TIER_PLACEHOLDER,
             cell_cycle_score_tbl_TIER_PLACEHOLDER,
             factors_to_regress = my_factors_to_regress
        ),
        pattern = map(slice(input_read, index  = SLICE_PLACEHOLDER ),
                      empty_droplets_tbl_TIER_PLACEHOLDER,
                      alive_identification_tbl_TIER_PLACEHOLDER,
                      cell_cycle_score_tbl_TIER_PLACEHOLDER),
        iteration = "list",
        resources = RESOURCE_PLACEHOLDER
      )}, 
      tiers = input_hpc$initialisation$tier,
      script = glue("{input_hpc$initialisation$store}.R")
    )
    
  }
  
  input_hpc |>
    c(list(normalise_abundance_seurat_SCT = args_list)) |>
    add_class("HPCell")
  
}

target_chunk_undefined_normalise_abundance_seurat_SCT = function(input_hpc){
  append_chunk_tiers(
    { tar_target(
      non_batch_variation_removal_S_TIER_PLACEHOLDER, 
      NULL, 
      pattern = slice(input_read, index  = SLICE_PLACEHOLDER ), 
      iteration = "list", 
      resources = RESOURCE_PLACEHOLDER
    ) }, 
    tiers = input_hpc$initialisation$tier,
    script = glue("{input_hpc$initialisation$store}.R")
  )
}

# Define the generic function
#' @export
calculate_pseudobulk <- function(input_data, group_by = NULL) {
  UseMethod("calculate_pseudobulk")
}

#' @export
calculate_pseudobulk.HPCell = function(input_hpc, group_by = NULL) {
  
  group_by = group_by |> enquo()
  
  # Capture all arguments including defaults
  args_list <- as.list(environment())[-1]
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  group_by |> saveRDS("pseudobulk_group_by.rds")
  
  args_list$target_chunk = function(){
    
    append_chunk_fix(
      { tar_target(pseudobulk_group_by, readRDS("pseudobulk_group_by.rds")) }, 
      script = glue("{input_hpc$initialisation$store}.R")
    )
    
    append_chunk_tiers(
      { tar_target(
        create_pseudobulk_sample_TIER_PLACEHOLDER, 
        create_pseudobulk(
          preprocessing_output_S_TIER_PLACEHOLDER,  
          sample_names,
          x = pseudobulk_group_by
        ), 
        pattern = map(preprocessing_output_S_TIER_PLACEHOLDER, slice(sample_names, index  = SLICE_PLACEHOLDER )), 
        iteration = "list",
        resources = RESOURCE_PLACEHOLDER
      )}, 
      tiers = input_hpc$initialisation$tier,
      script = glue("{input_hpc$initialisation$store}.R")
    )
    
  }

  input_hpc |>
    c(list(calculate_pseudobulk = args_list)) |>
    add_class("HPCell")
  
}

# Define the generic function
#' @export
evaluate_hpc <- function(input_data) {
  UseMethod("evaluate_hpc")
}

#' @importFrom glue glue
#' @importFrom purrr imap
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
  
    data_file_names |> as.list() |>  saveRDS("input_file.rds")
  
  input_hpc$initialisation$computing_resources |> saveRDS("temp_computing_resources.rds")
  tiers = input_hpc$initialisation$tier |> 
    get_positions() 
  tiers |> 
    saveRDS("temp_tiers.rds")
  
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
      controller = crew_controller_group ( readRDS("temp_computing_resources.rds") ), 
      packages = c("HPCell")
    )
    
    target_list = list(
      tar_target(read_file, readRDS("input_file.rds"), iteration = "list")
    )
    
    target_list = 
      target_list |> c(list(
        
        # Reading input files
        tar_target(
          input_read,
          readRDS(read_file),
          pattern = map(read_file),
          iteration = "list"
        ) ,
        tar_target(tiers, readRDS("temp_tiers.rds"))
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
    input_hpc$remove_empty_DropletUtils$target_chunk() 
    
  else
    target_chunk_undefined_remove_empty_DropletUtils(input_hpc)
  
  #-----------------------#
  # Annotate cell type
  #-----------------------#
  
  if(
    "annotate_cell_types" %in% names(input_hpc) |
    ( "remove_dead_scuttle" %in% names(input_hpc) & !is.null(input_hpc$remove_dead_scuttle$group_by))
  )
      input_hpc$annotate_cell_types$target_chunk()
    
  else
      target_chunk_undefined_annotate_cell_type(input_hpc)
  
  #-----------------------#
  # Remove dead
  #-----------------------#
  
  if("remove_dead_scuttle" %in% names(input_hpc))
    input_hpc$remove_dead_scuttle$target_chunk()

  else
      target_chunk_undefined_remove_dead_scuttle(input_hpc)
  
  
  #-----------------------#
  # score cell cycle
  #-----------------------#
  if("score_cell_cycle_seurat" %in% names(input_hpc))
    input_hpc$score_cell_cycle_seurat$target_chunk()

  else
      target_chunk_undefined_score_cell_cycle_seurat(input_hpc)
  
  #-----------------------#
  # Doublets
  #-----------------------#
  
  if("remove_doublets_scDblFinder" %in% names(input_hpc))
      input_hpc$remove_doublets_scDblFinder$target_chunk()
    
  else
    target_chunk_undefined_remove_doublets_scDblFinder(input_hpc)
  
  #-----------------------#
  # SCT
  #-----------------------#
  
  if("normalise_abundance_seurat_SCT" %in% names(input_hpc))
    input_hpc$normalise_abundance_seurat_SCT$target_chunk()
  
  else
    target_chunk_undefined_normalise_abundance_seurat_SCT(input_hpc)
  

  #-----------------------#
  # Create single cell output
  #-----------------------#
  
  append_chunk_tiers(
    { tar_target(
      preprocessing_output_S_TIER_PLACEHOLDER, 
      preprocessing_output(input_read,
           empty_droplets_tbl_TIER_PLACEHOLDER,
           non_batch_variation_removal_S_TIER_PLACEHOLDER,
           alive_identification_tbl_TIER_PLACEHOLDER,
           cell_cycle_score_tbl_TIER_PLACEHOLDER,
           annotation_label_transfer_tbl_TIER_PLACEHOLDER,
           doublet_identification_tbl_TIER_PLACEHOLDER),
      pattern = map(slice(input_read, index  = SLICE_PLACEHOLDER ),
                    empty_droplets_tbl_TIER_PLACEHOLDER,
                    non_batch_variation_removal_S_TIER_PLACEHOLDER,
                    alive_identification_tbl_TIER_PLACEHOLDER,
                    cell_cycle_score_tbl_TIER_PLACEHOLDER,
                    annotation_label_transfer_tbl_TIER_PLACEHOLDER,
                    doublet_identification_tbl_TIER_PLACEHOLDER),
      iteration = "list",
      resources = RESOURCE_PLACEHOLDER
    ) }, 
    tiers = input_hpc$initialisation$tier,
    script = glue("{input_hpc$initialisation$store}.R")
  )
  
  #-----------------------#
  # Create pseudobulk 
  #-----------------------#
  
  if("calculate_pseudobulk" %in% names(input_hpc))
      input_hpc$calculate_pseudobulk$target_chunk()

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
    tar_meta(store = glue("{input_hpc$initialisation$store}")) |> 
      filter(name |> str_detect("preprocessing_output_S_?[1-9]*$")) 
  )
}

#' @importFrom methods show
#' @export
print.HPCell <- function(x, ...){
  
  x |>
    evaluate_hpc() |> 
    print()
}
