#' Run Targets Pipeline for HPCell
#'
#' @description
#' This function sets up and executes a `targets` pipeline for HPCell. It saves input data and configurations,
#' writes a pipeline script, and runs the pipeline using the 'targets' package.
#'
#' @param input_hpc Character vector of input data path for the pipeline.
#' @param store Directory path for storing the pipeline files.
#' @param input_reference Optional reference data.
#' @param tissue Tissue type for the analysis.
#' @param computing_resources Configuration for computing resources.
#' @param debug_step Optional step for debugging.
#' @param filter_empty_droplets Flag to indicate if input filtering is needed.
#' @param RNA_assay_name Name of the RNA assay.
#' @param sample_column Column name for sample identification.
#' @param cell_type_annotation_column Column name for cell type annotation in input data
#' @param gene_nomenclature Character vector indicating gene nomenclature in input_data
#' @param data_container_type A character vector of length one specifies the input data type.
#' The accepted input data type are: 
#' sce_rds for `SingleCellExperiment` RDS,
#' seurat_rds for `Seurat` RDS,
#' sce_hdf5 for `SingleCellExperiment` HDF5-based object
#' seurat_h5 for `Seurat` HDF5-based object
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
initialise_hpc <- function(input_hpc,
                           store =  targets::tar_config_get("store"),
                           computing_resources = crew_controller_local(workers = 1),
                           tier = rep(1, length(input_hpc)),
                           debug_step = NULL,
                           RNA_assay_name = "RNA",
                           gene_nomenclature = "symbol",
                           data_container_type) {
  
  # Capture all arguments including defaults
  args_list <- as.list(environment())
  
  # if simple names are not set, use integers
  if(input_hpc |> names() |> is.null())
    input_hpc = input_hpc |> set_names(seq_len(length(input_hpc)))
  
  input_hpc |> names() |> saveRDS("sample_names.rds")
  #cell_count |> saveRDS("cell_count.rds")
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  # Write targets
  dir.create(store, showWarnings = FALSE, recursive = TRUE)
  data_file_names = glue("{store}/{names(input_hpc)}.rds")
  
  input_hpc |> as.list() |>  saveRDS("input_file.rds")
  gene_nomenclature |> saveRDS("temp_gene_nomenclature.rds")
  data_container_type |> saveRDS("data_container_type.rds")
  computing_resources |> saveRDS("temp_computing_resources.rds")
  tiers = tier |> 
    get_positions() 
  tiers |> 
    saveRDS("temp_tiers.rds")
  
  # Write pipeline to a file
  {
    library(HPCell)
    library(dplyr)
    library(magrittr)
    library(tibble)
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
      debug = d, # Set the target you want to debug.
      # cue = tar_cue(mode = "never") # Force skip non-debugging outdated targets.
      controller = crew_controller_group ( readRDS("temp_computing_resources.rds") ), 
      packages = c("HPCell", "tidySingleCellExperiment")
    )
    
    target_list = list(
      tar_target(gene_nomenclature, readRDS("temp_gene_nomenclature.rds"), iteration = "list", deployment = "main"),
      tar_target(data_container_type, readRDS("data_container_type.rds"), deployment = "main")
      
    )
    
    
    } |> 
    substitute(env = list(d = debug_step)) |> 
    tar_script_append2(script = glue("{store}.R"), append = FALSE)

  
  input_hpc = 
    list(initialisation = args_list ) |>
    c(list(sample_names = list(iterate = "map")) ) |> 
    
    add_class("HPCell")
  
  
  input_hpc |> 
    
    hpc_single(
      target_output = "sample_names", 
      user_function = readRDS |> quote(),
      file = "sample_names.rds", 
      deployment = "main", 
      iterate = "map"
    ) |> 
    
    hpc_single(
      target_output = "read_file_list", 
      user_function = readRDS |> quote(),
      file = "input_file.rds", 
      deployment = "main", 
      iterate = "map"
    ) |> 
  
    hpc_iterate(
      target_output = "data_object", 
      user_function = read_data_container |> quote() ,
      file  = "read_file_list" |> is_target(),
      container_type = data_container_type
    )
  
}





# Define the generic function
#' @export
remove_empty_DropletUtils <- function(input_hpc, total_RNA_count_check = NULL, target_input = "data_object", target_output = "empty_tbl", ...) {
  UseMethod("remove_empty_DropletUtils")
}

#' @export
remove_empty_DropletUtils.Seurat = function(input_hpc, total_RNA_count_check = NULL, target_input = "data_object", target_output = "empty_tbl", ...) {
  # Capture all arguments including defaults
  args_list <- as.list(environment())
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  list(initialisation = list(input_hpc = input_hpc)) |>
    add_class("HPCell") |>
    remove_empty_DropletUtils()
  
}

#' @export
remove_empty_DropletUtils.HPCell = function(input_hpc, total_RNA_count_check = NULL, target_input = "data_object", target_output = "empty_tbl",...) {

  input_hpc |> 
    hpc_iterate(
      target_output = target_output, 
      user_function = empty_droplet_id |> quote() , 
      input_read_RNA_assay = target_input |> is_target(),
      total_RNA_count_check = total_RNA_count_check
    )
  
}

target_chunk_undefined_remove_empty_DropletUtils = function(input_hpc){
  
  input_hpc |> 
    hpc_iterate(
      target_output = "empty_tbl", 
      user_function = (function(x) x |> as_tibble() |>  select(.cell) |> mutate(empty_droplet = FALSE)) |> quote() , 
      x = "data_object" |> is_target(),
      packages = c("dplyr", "tidySingleCellExperiment", "tidyseurat")
    )
  
}


# Define the generic function
#' @export
remove_dead_scuttle <- function(input_hpc, 
                                group_by = NULL, 
                                target_output = "alive_tbl",
                                target_input = "data_object", 
                                target_empty_droplets = "empty_tbl",
                                target_annotation = NULL) {
  UseMethod("remove_dead_scuttle")
}

#' @export

remove_dead_scuttle.HPCell = function(
    input_hpc, 
    group_by = NULL, 
    target_output = "alive_tbl",
    target_input = "data_object", 
    target_empty_droplets = "empty_tbl",
    target_annotation = NULL
    
  ) {
  
  input_hpc |> 
    hpc_iterate(
      target_output = target_output, 
      user_function = alive_identification |> quote() , 
      input_read_RNA_assay = target_input |> safe_as_name(), 
      empty_droplets_tbl = target_empty_droplets |> safe_as_name() ,
      annotation_label_transfer_tbl = target_annotation |> safe_as_name() ,
      annotation_column = group_by 
    )
  
}

#' @importFrom dplyr mutate
target_chunk_undefined_remove_dead_scuttle = function(input_hpc){
  
  input_hpc |> 
    hpc_iterate(
      target_output = "alive_tbl", 
      user_function = (function(x) x |> as_tibble() |>  select(.cell) |> mutate(alive = TRUE)) |> quote() , 
      x = "data_object" |> is_target(),
      packages = c("dplyr", "tidySingleCellExperiment", "tidyseurat")
    )
  
}



# Define the generic function
#' @export
score_cell_cycle_seurat <- function(input_hpc, target_input = "data_object", target_output = "cell_cycle_tbl",...) {
  UseMethod("score_cell_cycle_seurat")
}

#' @export
score_cell_cycle_seurat.HPCell = function(input_hpc, target_input = "data_object", target_output = "cell_cycle_tbl", ...) {
  
  input_hpc |> 
    hpc_iterate(
      target_output = target_output, 
      user_function = cell_cycle_scoring |> quote() , 
      input_read_RNA_assay = target_input |> is_target(), 
      empty_droplets_tbl = "empty_tbl" |> is_target() ,
      gene_nomenclature = input_hpc$initialisation$gene_nomenclature
    )
  
}

target_chunk_undefined_score_cell_cycle_seurat = function(input_hpc, target_input = "data_object"){
  
  input_hpc |> 
    hpc_iterate(
      target_output = "cell_cycle_tbl", 
      user_function = function(x) NULL  , 
      x = "read_file_list" |> is_target()
    )
  
}


# Define the generic function
#' @export
remove_doublets_scDblFinder <- function(input_hpc, target_input = "data_object", target_output = "doublet_tbl") {
  UseMethod("remove_doublets_scDblFinder")
}

#' @export
remove_doublets_scDblFinder.HPCell = function(input_hpc, target_input = "data_object", target_output = "doublet_tbl") {
  input_hpc |> 
    hpc_iterate(
      target_output = target_output, 
      user_function = doublet_identification |> quote() , 
      input_read_RNA_assay = target_input |> is_target(), 
      empty_droplets_tbl = "empty_tbl" |> is_target() ,
      alive_identification_tbl = "alive_tbl" |> is_target()
    )
  
}

target_chunk_undefined_remove_doublets_scDblFinder = function(input_hpc){
  
  input_hpc |> 
    hpc_iterate(
      target_output = "doublet_tbl", 
      user_function = (function(x) x |> as_tibble() |>  select(.cell) |> mutate(scDblFinder.class="singlet")) |> quote() , 
      x = "data_object" |> is_target(),
      packages = c("dplyr", "tidySingleCellExperiment", "tidyseurat")
    )
  
}

# Define the generic function
#' @export
annotate_cell_type <- function(input_hpc, azimuth_reference = NULL, target_input = "data_object", target_output = "annotation_tbl",...) {
  UseMethod("annotate_cell_type")
}

#' @export
annotate_cell_type.HPCell = function(input_hpc, azimuth_reference = NULL, target_input = "data_object", target_output = "annotation_tbl", ...) {
  
  azimuth_reference |> saveRDS("input_reference.rds")
  
  input_hpc |> 
    
    hpc_single(
      target_output = "reference_read", 
      user_function = readRDS |> quote(),
      file = "input_reference.rds" 
    ) |> 
  
    hpc_iterate(
      target_output = target_output, 
      user_function = annotation_label_transfer |> quote() , 
      input_read_RNA_assay = target_input |> is_target(), 
      empty_droplets_tbl = "empty_tbl" |> is_target() ,
      reference_azimuth = reference_read |> quote()
    )
  
}

target_chunk_undefined_annotate_cell_type = function(input_hpc){
  
  input_hpc |> 
    hpc_iterate(
      target_output = "annotation_tbl", 
      user_function = function(x) NULL ,
      x = read_file_list |> quote()
    )
  
 
}


# Define the generic function
#' @export
normalise_abundance_seurat_SCT <- function(input_hpc, target_input = "data_object", target_output = "sct_matrix", ...) {
  UseMethod("normalise_abundance_seurat_SCT")
}

#' @export
normalise_abundance_seurat_SCT.HPCell = function(input_hpc, factors_to_regress = NULL, target_input = "data_object", target_output = "sct_matrix", ...) {
  
  input_hpc |> 
  hpc_iterate(
    target_output = target_output, 
    user_function = non_batch_variation_removal |> quote() , 
    input_read_RNA_assay = target_input |> is_target(), 
    empty_droplets_tbl = "empty_tbl" |> is_target() ,
    alive_identification_tbl = "alive_tbl" |> is_target(),
    cell_cycle_score_tbl = "cell_cycle_tbl" |> is_target(),
    factors_to_regress = factors_to_regress,
    external_path = glue("{input_hpc$initialisation$store}/external"),
    ...
  )
  
}

target_chunk_undefined_normalise_abundance_seurat_SCT = function(input_hpc){
  
  input_hpc |> 
    hpc_iterate(
      target_output = "sct_matrix", 
      user_function = function(x) NULL ,
      x = read_file_list |> quote()
    )
  

}

# Define the generic function
#' @export
calculate_pseudobulk <- function(input_hpc, group_by = NULL, target_input = "data_object", target_output = "pseudobulk_se") {
  UseMethod("calculate_pseudobulk")
}

#' @export
calculate_pseudobulk.HPCell = function(input_hpc, group_by = NULL, target_input = "data_object", target_output = "pseudobulk_se") {
  
  pseudobulk_sample = 
    glue("{target_output}_iterated") |> 
    # This is important otherwise targets fails with glue
    as.character()
  
  input_hpc |> 
    
    # aggregate each
    hpc_iterate(
      target_output = pseudobulk_sample, 
      user_function = create_pseudobulk |> quote() , 
      input_read_RNA_assay = target_input |> is_target(), 
      sample_names_vec = "sample_names" |> is_target(),
      empty_droplets_tbl = "empty_tbl" |> is_target() ,
      alive_identification_tbl = "alive_tbl" |> is_target(),
      cell_cycle_score_tbl = "cell_cycle_tbl" |> is_target(),
      annotation_label_transfer_tbl = "annotation_tbl" |> is_target(),
      doublet_identification_tbl = "doublet_tbl" |> is_target(),
      x = group_by,
      external_path = glue("{input_hpc$initialisation$store}/external")
    ) |> 
    
    # merge
    hpc_merge(
      target_output = target_output, 
      user_function = pseudobulk_merge |> quote(), 
      external_path = glue("{input_hpc$initialisation$store}/external"),
      pseudobulk_list = pseudobulk_sample |> is_target(),
      packages = c("tidySummarizedExperiment", "HPCell")
    )
  
  
}

# Define the generic function
#' @export
get_single_cell <- function(input_hpc, target_input = "data_object", target_output = "single_cell",...) {
  UseMethod("get_single_cell")
}

#' @export
get_single_cell.HPCell = function(input_hpc, target_input = "data_object", target_output = "single_cell", ...) {
  
  input_hpc |> 
    
    # aggregate each
    hpc_iterate(
      target_output = target_output, 
      user_function = preprocessing_output |> quote() , 
      input_read_RNA_assay = target_input |> is_target(), 
      empty_droplets_tbl = "empty_tbl" |> is_target() ,
      non_batch_variation_removal_S = sct_matrix |> quote(), 
      alive_identification_tbl = "alive_tbl" |> is_target(),
      cell_cycle_score_tbl = "cell_cycle_tbl" |> is_target(),
      annotation_label_transfer_tbl = "annotation_tbl" |> is_target(),
      doublet_identification_tbl = "doublet_tbl" |> is_target()
    )
  
  
}


#' Test Differential Abundance for HPCell
#'
#' This function tests differential abundance for HPCell objects.
#'
#' @name test_differential_abundance,HPCell-method
#' @rdname test_differential_abundance
#' @inherit tidybulk::test_differential_abundance
#'
#' @importFrom tidybulk test_differential_abundance
#' @exportMethod test_differential_abundance
#' @param .data An HPCell object.
#' @param .formula A formula used to model the design matrix.
#' @param .sample Sample parameter.
#' @param .transcript Transcript parameter.
#' @param .abundance Abundance parameter.
#' @param contrasts Contrasts parameter.
#' @param method Method parameter, default is "edgeR_quasi_likelihood".
#' @param test_above_log2_fold_change Test above log2 fold change.
#' @param scaling_method Scaling method, default is "TMM".
#' @param omit_contrast_in_colnames Omit contrast in column names.
#' @param prefix Prefix parameter.
#' @param action Action parameter, default is "add".
#' @param ... Additional parameters.
#' @param significance_threshold Significance threshold.
#' @param fill_missing_values Fill missing values.
#' @param .contrasts Contrasts parameter.
#' @return The result of the differential abundance test.
#'
setMethod(
  "test_differential_abundance",
  signature(.data = "HPCell"),
  function(.data, .formula, .sample = NULL, .transcript = NULL, 
           .abundance = NULL, contrasts = NULL, method = "edgeR_quasi_likelihood", 
           test_above_log2_fold_change = NULL, scaling_method = "TMM", 
           omit_contrast_in_colnames = FALSE, prefix = "", action = "add", factor_of_interest = NULL,
           target_input = "pseudobulk_se", target_output = "de", group_by_column = NULL,
           ..., significance_threshold = NULL, fill_missing_values = NULL, 
           .contrasts = NULL) {
    
      if(.formula |> deparse() |> str_detect("\\|"))
        factory_de_random_effect(
          se_list_input = target_input, 
          output_se = target_output, 
          formula=.formula,
          #method="edger_robust_likelihood_ratio", 
          tiers = tiers,
          factor_of_interest = factor_of_interest,
          .abundance = .abundance
        )
      
      else
       
        .data |> 
          
          hpc_single(
            target_output = "chunk_tbl", 
            user_function = function(x){ x |> rownames() |> feature_chunks()} |> quote(),
            x = "pseudobulk_se" |> is_target() 
          ) |> 
        
          hpc_single(
            target_output = "pseudobulk_group_list", 
            user_function = group_split |> quote(),
            .tbl  = target_input  |> is_target(),
            gr = as.name(gr) |> substitute(env = list(gr = group_by_column)),
            packages = c("tidySummarizedExperiment", "S4Vectors", "targets"), 
            
            # I need this because targets does not know the output 
            # is a list I need to iterate on outside the tiers
            iterate = "map"
          ) |> 
          

          hpc_iterate(
            target_output = target_output, 
            user_function = internal_de_function |> quote() , 
            x = "pseudobulk_group_list" |> is_target(),
            fi = factor_of_interest,
            a = .abundance,
            f = .formula,
            m = method,
            packages="tidybulk"
          ) 
        
      
        
})



    #   
    #   factory_collapse(
    #     "colapsed_preprocessing_output",
    #     bind_rows(preprocessing_output_S) ,
    #     "preprocessing_output_S",
    #     tiers
    #   ),
    #   
    #   tar_render(
    #     name = preprocessing_report,
    #     path = paste0(system.file(package = "HPCell"), "/rmd/preprocessing_report.Rmd"),
    #     params = list(
    #       x1 = collapsed_preprocessing_output, 
    #       x2 = group_by|> quo_name()
    #     )
    # ) 


#' generate_report = function(tiers){
#' 
#'   list(
#'     factory_split(
#'       "final_report",
#'       command = {read_file |>
#'         read_data_container(container_type = data_container_type) |>
        # tar_render(
        #   name = empty_droplets_report,
        #   path =  paste0(system.file(package = "HPCell"), "/rmd/Empty_droplet_report.Rmd"),
        #   params = list(x1 = empty_droplets_tbl,
        #                 # x2 = empty_droplets_tbl,
        #                 # x3 = annotation_label_transfer_tbl
        #                 # x4 = tar_read(unique_tissues, store = store),
        #                 # x5 = sample_column |> quo_name()
        #                 )) |>
#'           quote()
#'         },
#'         tiers,
#'         arguments_to_tier = "read_file",
#'         other_arguments_to_tier = c("empty_droplets_tbl"
#'                                     # "annotation_label_transfer_tbl",
#'                                     # "doublet_identification_tbl"),
#'         ),
#'         other_arguments_to_map = c("empty_droplets_tbl"
#'                                    # "annotation_label_transfer_tbl",
#'                                    # "doublet_identification_tbl")
#'       )
#'     )
#' 
#'   )
#' 
#' }



# Define the generic function
#' @export
evaluate_hpc <- function(input_hpc) {
  UseMethod("evaluate_hpc")
}

#' @importFrom glue glue
#' @importFrom purrr imap
#' @export
evaluate_hpc.HPCell = function(input_hpc) {
  
  #-----------------------#
  # Empty droplets
  #-----------------------#
  
  if(! "empty_tbl" %in% names(input_hpc))
    target_chunk_undefined_remove_empty_DropletUtils(input_hpc)
  
  #-----------------------#
  # Annotate cell type
  #-----------------------#
  
  if(
    !("annotation_tbl" %in% names(input_hpc) |
      ( "alive_tbl" %in% names(input_hpc) & !is.null(input_hpc$remove_dead_scuttle$group_by))
    ))
    target_chunk_undefined_annotate_cell_type(input_hpc)
  
  #-----------------------#
  # Remove dead
  #-----------------------#
  
  if(! "alive_tbl" %in% names(input_hpc))
    target_chunk_undefined_remove_dead_scuttle(input_hpc)
  
  
  #-----------------------#
  # score cell cycle
  #-----------------------#
  if(! "cell_cycle_tbl" %in% names(input_hpc))
    target_chunk_undefined_score_cell_cycle_seurat(input_hpc)
  
  #-----------------------#
  # Doublets
  #-----------------------#
  if(! "doublet_tbl" %in% names(input_hpc))
    target_chunk_undefined_remove_doublets_scDblFinder(input_hpc)
  
  #-----------------------#
  # SCT
  #-----------------------#
  
  if(! "sct_matrix" %in% names(input_hpc))
    target_chunk_undefined_normalise_abundance_seurat_SCT(input_hpc)
  
  
  #-----------------------#
  # Reports
  #-----------------------#
  # tar_tier_append(
  #   quote(generate_report),
  #   input_hpc$initialisation$tier |> get_positions() ,
  #   glue("{input_hpc$initialisation$store}.R")
  # )
  
  #-----------------------#
  # Close pipeline
  #-----------------------#
  
  # Call final list
  tar_script_append({
    target_list
  }, script = glue("{input_hpc$initialisation$store}.R"))
  
  if(input_hpc$initialisation$debug_step |> is.null())
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
    "pseudobulk_group_by.rds",
    "temp_tiers.rds",
    "temp_gene_nomenclature.rds",
    "data_container_type.rds"
  ) |> 
    remove_files_safely()
  
  # If get_single_cell is called then return the object 
  if(input_hpc$last_call |> is.null() |> not())
    return(
      tar_meta(store = glue("{input_hpc$initialisation$store}")) |> 
        filter(name |> str_detect("single_cell_?.*$"), type=="pattern") |> 
        pull(name) |> 
        map(tar_read_raw) |> 
        unlist() |> 
        do.call(cbind, args = _) 
    )
  
  else
    return(
      tar_meta(store = glue("{input_hpc$initialisation$store}")) |> 
        filter(name |> str_detect("single_cell_?.*$"), type=="pattern") 
    )
}

#' @importFrom methods show
#' @export
print.HPCell <- function(x, ...){
  
  x |>
    evaluate_hpc() |> 
    print()
}