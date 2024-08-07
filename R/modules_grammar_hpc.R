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
  # map2(
  #   input_hpc,
  #   data_file_names,
  #   ~ .x |> saveRDS(.y)
  # )
  
  input_hpc |> as.list() |>  saveRDS("input_file.rds")
  gene_nomenclature |> saveRDS("temp_gene_nomenclature.rds")
  data_container_type |> saveRDS("data_container_type.rds")
  
  computing_resources |> saveRDS("temp_computing_resources.rds")
  tiers = tier |> 
    get_positions() 
  tiers |> 
    saveRDS("temp_tiers.rds")
  
  #sample_column |> saveRDS("sample_column.rds")
  debug_step |> saveRDS("temp_debug_step.rds")
  
  # Write pipeline to a file
  tar_script({
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
      debug = readRDS("temp_debug_step.rds"), # Set the target you want to debug.
      # cue = tar_cue(mode = "never") # Force skip non-debugging outdated targets.
      controller = crew_controller_group ( readRDS("temp_computing_resources.rds") ), 
      packages = c("HPCell")
    )
    
    target_list = list(
      tar_files(read_file,readRDS("input_file.rds") , iteration = "list", deployment = "main"),
      #tar_target(read_file, readRDS("input_file.rds"), format = "file", iteration = "list"),
      tar_target(gene_nomenclature, readRDS("temp_gene_nomenclature.rds"), iteration = "list", deployment = "main"),
      tar_target(data_container_type, readRDS("data_container_type.rds"), deployment = "main")
    
    )
    
    target_list = 
      target_list |> c(list(
        
        # # Reading input files
        # tar_target(
        #   input_read,
        #   readRDS(read_file),
        #   pattern = map(read_file),
        #   iteration = "list"
        # ) ,
        tar_target(tiers, readRDS("temp_tiers.rds"), deployment = "main")
      ))
    
  }, script = glue("{store}.R"), ask = FALSE)
  
  # Set sample names
  tar_script_append2(
    substitute({
      
      target_list = 
        target_list |> c(list(
          tar_target(sample_names, readRDS("sample_names.rds"))
        ))
      
    }),
    script = glue("{store}.R")
  )
  
  
  list(initialisation = args_list) |>
    add_class("HPCell")
}





# Define the generic function
#' @export
remove_empty_DropletUtils <- function(input_hpc, total_RNA_count_check = NULL, ...) {
  UseMethod("remove_empty_DropletUtils")
}

#' @export
remove_empty_DropletUtils.Seurat = function(input_hpc, total_RNA_count_check = NULL, ...) {
  # Capture all arguments including defaults
  args_list <- as.list(environment())
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  list(initialisation = list(input_hpc = input_hpc)) |>
    add_class("HPCell") |>
    remove_empty_DropletUtils()
  
}

#' @export
remove_empty_DropletUtils.HPCell = function(input_hpc, total_RNA_count_check = NULL, ...) {
  
  # Capture all arguments including defaults
  args_list <- as.list(environment())[-1]
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  args_list$factory = function(tiers){
    list(
      tar_target_raw("total_RNA_count_check", readRDS("total_RNA_count_check.rds") |> quote(), deployment = "main"),
      
      factory_split(
        "empty_droplets_tbl", 
        read_file |> 
          read_data_container(container_type = data_container_type) |> 
          empty_droplet_id(total_RNA_count_check) |> 
          quote(),
        tiers, 
        arguments_to_tier = "read_file"
      )  ,
      
      factory_collapse(
        "my_report",
        bind_rows(empty_droplets_tbl) |> quote(),
        "empty_droplets_tbl",
        tiers, packages = c("dplyr")
      )
    )
    
  }
  
  # We don't want recursive when we call factory
  if(input_hpc |> length() > 0) {
    total_RNA_count_check |> saveRDS("total_RNA_count_check.rds")
    
    tar_tier_append(
      quote(dummy_hpc |> remove_empty_DropletUtils() %$% remove_empty_DropletUtils %$% factory),
      input_hpc$initialisation$tier |> get_positions() ,
      glue("{input_hpc$initialisation$store}.R")
    )
    
  }

  
  # Add pipeline step
  input_hpc |>
    c(list(remove_empty_DropletUtils = args_list)) |>
    add_class("HPCell")
  
  
}

target_chunk_undefined_remove_empty_DropletUtils = function(input_hpc){
  append_chunk_tiers(
    { tar_target(
      empty_droplets_tbl_TIER_PLACEHOLDER,
      read_file |> 
        read_data_container(container_type = data_container_type) |> as_tibble() |> select(.cell) |> mutate(empty_droplet = FALSE),
      pattern = slice(read_file, index  = SLICE_PLACEHOLDER ),
      iteration = "list", 
      resources = RESOURCE_PLACEHOLDER,
      packages = c("dplyr", "tidySingleCellExperiment", "tidyseurat")
    ) }, 
    tiers = input_hpc$initialisation$tier,
    script = glue("{input_hpc$initialisation$store}.R")
  )
}


# Define the generic function
#' @export
remove_dead_scuttle <- function(input_hpc, group_by = NULL) {
  UseMethod("remove_dead_scuttle")
}

#' @export
remove_dead_scuttle.HPCell = function(input_hpc, group_by = NULL) {
  
  # Capture all arguments including defaults
  args_list <- as.list(environment())[-1]
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  args_list$factory = function(tiers){
    list(
      tar_target_raw("grouping_column", readRDS("temp_group_by.rds") |> quote(), deployment = "main"),
      
      factory_split(
        "alive_identification_tbl", 
        read_file |> 
          read_data_container(container_type = data_container_type) |> 
          alive_identification(
            empty_droplets_tbl,
            annotation_label_transfer_tbl,
            grouping_column
          ) |> quote(),
        tiers, arguments_to_tier = "read_file",
        other_arguments_to_tier = c("empty_droplets_tbl", "annotation_label_transfer_tbl"), 
        other_arguments_to_map = c("empty_droplets_tbl", "annotation_label_transfer_tbl")
      ),
      
      factory_collapse(
        "my_report2",
        bind_rows(alive_identification_tbl) |> quote(), 
        "alive_identification_tbl",
        tiers, packages = c("dplyr")
      )
    )
    
  }
  
  # We don't want recursive when we call factory
  if(input_hpc |> length() > 0) {
    group_by |> saveRDS("temp_group_by.rds")
    
    tar_tier_append(
      quote(dummy_hpc |> remove_dead_scuttle() %$% remove_dead_scuttle %$% factory),
      input_hpc$initialisation$tier |> get_positions() ,
      glue("{input_hpc$initialisation$store}.R")
    )
  }

  
  input_hpc |>
    c(list(remove_dead_scuttle = args_list))  |>
    add_class("HPCell")
  
}

#' @importFrom dplyr mutate
target_chunk_undefined_remove_dead_scuttle = function(input_hpc){
  append_chunk_tiers(
    { tar_target(
      alive_identification_tbl_TIER_PLACEHOLDER, 
      read_file |> 
        read_data_container(container_type = data_container_type) |> as_tibble() |> dplyr::select(.cell) |> mutate(alive = TRUE), 
      pattern = slice(read_file, index  = SLICE_PLACEHOLDER ), 
      iteration = "list", 
      resources = RESOURCE_PLACEHOLDER,
      packages = c("dplyr", "tidySingleCellExperiment", "tidyseurat")
    ) }, 
    tiers = input_hpc$initialisation$tier,
    script = glue("{input_hpc$initialisation$store}.R")
  )
}



# Define the generic function
#' @export
score_cell_cycle_seurat <- function(input_hpc, ...) {
  UseMethod("score_cell_cycle_seurat")
}

#' @export
score_cell_cycle_seurat.HPCell = function(input_hpc) {
  # Capture all arguments including defaults
  args_list <- as.list(environment())[-1]
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  # Add pipeline step
  args_list$factory = function(tiers){
    list(
      factory_split(
        "cell_cycle_score_tbl", read_file |> 
          read_data_container(container_type = data_container_type) |> 
          cell_cycle_scoring(
            empty_droplets_tbl,  
            gene_nomenclature = gene_nomenclature
          ) |> quote(),
        tiers, arguments_to_tier = "read_file",
        other_arguments_to_tier = c("empty_droplets_tbl"), other_arguments_to_map = c("empty_droplets_tbl")
      ),
      
      factory_collapse(
        "my_report3",
        bind_rows(cell_cycle_score_tbl) |> quote(), 
        "cell_cycle_score_tbl",
        tiers, packages = c("dplyr")
      )
    )
    
  }
  
  # We don't want recursive when we call factory
  if(input_hpc |> length() > 0) 
    tar_tier_append(
      quote(dummy_hpc |> score_cell_cycle_seurat() %$% score_cell_cycle_seurat %$% factory),
      input_hpc$initialisation$tier |> get_positions() ,
      glue("{input_hpc$initialisation$store}.R")
    )
  
  input_hpc |>
    c(list(score_cell_cycle_seurat = args_list)) |>
    add_class("HPCell")
  
}

target_chunk_undefined_score_cell_cycle_seurat = function(input_hpc){
  append_chunk_tiers(
    { tar_target(
      cell_cycle_score_tbl_TIER_PLACEHOLDER,
      NULL,
      pattern = slice(read_file, index  = SLICE_PLACEHOLDER ),
      iteration = "list",
      resources = RESOURCE_PLACEHOLDER,
      packages = c("dplyr", "tidySingleCellExperiment", "tidyseurat")
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
  
  args_list$factory = function(tiers){
    list(
      factory_split(
        "doublet_identification_tbl", read_file |> 
          read_data_container(container_type = data_container_type) |> 
          doublet_identification(
            empty_droplets_tbl,
            alive_identification_tbl
          ) |> quote(),
        tiers, arguments_to_tier = "read_file",
        other_arguments_to_tier = c("empty_droplets_tbl", "alive_identification_tbl"), other_arguments_to_map = c("empty_droplets_tbl", "alive_identification_tbl")
      ),
      
      factory_collapse(
        "my_report4",
        bind_rows(doublet_identification_tbl) |> quote(), 
        "doublet_identification_tbl",
        tiers, packages = c("dplyr")
      )
    )
    
  }
  
  # We don't want recursive when we call factory
  if(input_hpc |> length() > 0) 
    tar_tier_append(
      quote(dummy_hpc |> remove_doublets_scDblFinder() %$% remove_doublets_scDblFinder %$% factory),
      input_hpc$initialisation$tier |> get_positions() ,
      glue("{input_hpc$initialisation$store}.R")
    )
  
  input_hpc |>
    c(list(remove_doublets_scDblFinder = args_list)) |>
    add_class("HPCell")
  
}

target_chunk_undefined_remove_doublets_scDblFinder = function(input_hpc){
  append_chunk_tiers(
    { tar_target(
      doublet_identification_tbl_TIER_PLACEHOLDER, 
      read_file |> 
        read_data_container(container_type = data_container_type) |> as_tibble() |> select(.cell) |> mutate(scDblFinder.class="singlet"), 
      pattern = slice(read_file, index  = SLICE_PLACEHOLDER ), 
      iteration = "list", 
      resources = RESOURCE_PLACEHOLDER,
      packages = c("dplyr", "tidySingleCellExperiment", "tidyseurat")
    ) }, 
    tiers = input_hpc$initialisation$tier,
    script = glue("{input_hpc$initialisation$store}.R")
  )
}

# Define the generic function
#' @export
annotate_cell_type <- function(input_hpc, ...) {
  UseMethod("annotate_cell_type")
}

#' @export
annotate_cell_type.HPCell = function(input_hpc, azimuth_reference = NULL) {
  # Capture all arguments including defaults
  args_list <- as.list(environment())[-1]
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  args_list$factory = function(tiers){
    list(
      tar_target_raw("reference_read", readRDS("input_reference.rds") |> quote()),
      
      factory_split(
        "annotation_label_transfer_tbl", read_file |> 
          read_data_container(container_type = data_container_type) |> 
          annotation_label_transfer(
            empty_droplets_tbl,
            reference_read
          ) |> quote(),
        tiers, arguments_to_tier = "read_file",
        other_arguments_to_tier = c("empty_droplets_tbl"), other_arguments_to_map = c("empty_droplets_tbl") 
      ),
      
      factory_collapse(
        "my_report5",
        bind_rows(annotation_label_transfer_tbl) |> quote(), 
        "annotation_label_transfer_tbl",
        tiers, packages = c("dplyr")
      )
    )
    
  }
  
  # We don't want recursive when we call factory
  if(input_hpc |> length() > 0) {
    azimuth_reference |> saveRDS("input_reference.rds")
    
    tar_tier_append(
      quote(dummy_hpc |> annotate_cell_type() %$% annotate_cell_type %$% factory),
      input_hpc$initialisation$tier |> get_positions() ,
      glue("{input_hpc$initialisation$store}.R")
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
      pattern = slice(read_file, index  = SLICE_PLACEHOLDER ), 
      iteration = "list", 
      resources = RESOURCE_PLACEHOLDER,
      packages = c("dplyr", "tidySingleCellExperiment", "tidyseurat")
    ) }, 
    tiers = input_hpc$initialisation$tier,
    script = glue("{input_hpc$initialisation$store}.R")
  )
}


# Define the generic function
#' @export
normalise_abundance_seurat_SCT <- function(input_hpc, ...) {
  UseMethod("normalise_abundance_seurat_SCT")
}

#' @export
normalise_abundance_seurat_SCT.HPCell = function(input_hpc, factors_to_regress = NULL) {
  # Capture all arguments including defaults
  args_list <- as.list(environment())[-1]
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  args_list$factory = function(tiers){
    list(
      tar_target_raw("factors_to_regress", readRDS("factors_to_regress.rds") |> quote(), deployment = "main"),
      
      factory_split(
        "non_batch_variation_removal_S", read_file |> 
          read_data_container(container_type = data_container_type) |> 
          non_batch_variation_removal(
            empty_droplets_tbl,
            alive_identification_tbl,
            cell_cycle_score_tbl,
            factors_to_regress = factors_to_regress
          ) |> quote(),
        tiers, arguments_to_tier = "read_file",
        other_arguments_to_tier = c("empty_droplets_tbl", "alive_identification_tbl", "cell_cycle_score_tbl"), other_arguments_to_map = c("empty_droplets_tbl", "alive_identification_tbl", "cell_cycle_score_tbl")
        
      )
      # ,
      # 
      # factory_collapse(
      #   "my_report6",
      #   bind_rows(tibble()) |> quote(), 
      #   "non_batch_variation_removal_S",
      #   tiers, packages = c("dplyr", "tidySingleCellExperiment", "tidyseurat")
      # )
    )
    
  }
  
  # We don't want recursive when we call factory
  if(input_hpc |> length() > 0) {
    factors_to_regress |> saveRDS("factors_to_regress.rds")
    
    tar_tier_append(
      quote(dummy_hpc |> normalise_abundance_seurat_SCT() %$% normalise_abundance_seurat_SCT %$% factory),
      input_hpc$initialisation$tier |> get_positions() ,
      glue("{input_hpc$initialisation$store}.R")
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
      pattern = slice(read_file, index  = SLICE_PLACEHOLDER ), 
      iteration = "list", 
      resources = RESOURCE_PLACEHOLDER,
      packages = c("dplyr", "tidySingleCellExperiment", "tidyseurat")
    ) }, 
    tiers = input_hpc$initialisation$tier,
    script = glue("{input_hpc$initialisation$store}.R")
  )
}

# Define the generic function
#' @export
calculate_pseudobulk <- function(input_hpc, group_by = NULL) {
  UseMethod("calculate_pseudobulk")
}

#' @export
calculate_pseudobulk.HPCell = function(input_hpc, group_by = NULL) {
  
  # Capture all arguments including defaults
  args_list <- as.list(environment())[-1]
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  args_list$factory = function(tiers, external_path, pseudobulk_group_by = ""){
    
    list(
      tar_target_raw("pseudobulk_group_by", pseudobulk_group_by, deployment = "main") ,

      factory_split(
        "create_pseudobulk_sample", 
        read_file |> 
          read_data_container(container_type = data_container_type) |> 
          create_pseudobulk(
            sample_names, 
            empty_droplets_tbl,
            alive_identification_tbl,
            cell_cycle_score_tbl,
            annotation_label_transfer_tbl,
            doublet_identification_tbl,
            
            x = pseudobulk_group_by, 
            external_path = e
          ) |> 
          substitute(env = list(e = external_path)),
        tiers, arguments_to_tier = c("read_file", "sample_names"), 
        other_arguments_to_tier = c("empty_droplets_tbl",
                                    "alive_identification_tbl",
                                    "cell_cycle_score_tbl",
                                    "annotation_label_transfer_tbl",
                                    "doublet_identification_tbl"),
        other_arguments_to_map = c("empty_droplets_tbl",
                                   "alive_identification_tbl",
                                   "cell_cycle_score_tbl",
                                   "annotation_label_transfer_tbl",
                                   "doublet_identification_tbl")
        
      ),
      factory_merge_pseudobulk(
        se_list_input = "create_pseudobulk_sample",
        "pseudobulk_gran_group", 
        tiers, 
        external_path = external_path
      ) 
      
      # ,
      # 
      # factory_collapse(
      #   "my_report7",
      #   bind_rows(create_pseudobulk_sample) |> quote(),
      #   "create_pseudobulk_sample",
      #   tiers
      # )
    )
    
  }
   
  # We don't want recursive when we call factory
  if(input_hpc |> length() > 0) {

    tar_tier_append(
      quote(dummy_hpc |> calculate_pseudobulk() %$% calculate_pseudobulk %$% factory),
      input_hpc$initialisation$tier |> get_positions() ,
      script = glue("{input_hpc$initialisation$store}.R"), 
      external_path = glue("{input_hpc$initialisation$store}/external"),
      pseudobulk_group_by = group_by
    )
  }
  
  
  input_hpc |>
    c(list(calculate_pseudobulk = args_list)) |>
    add_class("HPCell")
  
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
           ..., significance_threshold = NULL, fill_missing_values = NULL, 
           .contrasts = NULL) {
    
    # Capture all arguments including defaults
    args_list <- as.list(environment())[-1]
    
    # Optionally, you can evaluate the arguments if they are expressions
    args_list <- lapply(args_list, eval, envir = parent.frame())
    
    args_list$factory = function(tiers, .formula, factor_of_interest = NULL, .abundance = NULL){
        
      if(.formula |> deparse() |> str_detect("\\|"))
        factory_de_random_effect(
          se_list_input = "create_pseudobulk_sample", 
          output_se = "de", 
          formula=.formula,
          #method="edger_robust_likelihood_ratio", 
          tiers = tiers,
          factor_of_interest = factor_of_interest,
          .abundance = .abundance
        )
      
      else
        factory_de_fix_effect(
          se_list_input = "create_pseudobulk_sample", 
          output_se = "de", 
          formula=.formula,
          method="edger_robust_likelihood_ratio", 
          tiers = tiers,
          factor_of_interest = factor_of_interest,
          .abundance = .abundance
        )
      
    }
    
    # We don't want recursive when we call factory
    if(.data |> length() > 0) {
      
      environment(.formula) <- new.env(parent = emptyenv())
      
      tar_tier_append(
        quote(dummy_hpc |> test_differential_abundance() %$% test_differential_abundance %$% factory),
        tiers = .data$initialisation$tier |> get_positions() ,
        script = glue("{.data$initialisation$store}.R"),
        .formula = .formula, 
        factor_of_interest = factor_of_interest,
        .abundance = .abundance
      )
      
    }
    
    
    .data |>
      c(list(test_differential_abundance = args_list)) |>
      add_class("HPCell")
    
  }
)


#' @export
preprocessing_output_factory = function(tiers){
  
  
  
  list(
    factory_split(
      "preprocessing_output_S", 
      read_file |> 
        read_data_container(container_type = data_container_type) |> 
        preprocessing_output(empty_droplets_tbl,
                             non_batch_variation_removal_S,
                             alive_identification_tbl,
                             cell_cycle_score_tbl,
                             annotation_label_transfer_tbl,
                             doublet_identification_tbl) |> 
        quote(),
      tiers, 
      arguments_to_tier = "read_file", 
      other_arguments_to_tier = c("empty_droplets_tbl",
                                  "non_batch_variation_removal_S",
                                  "alive_identification_tbl",
                                  "cell_cycle_score_tbl",
                                  "annotation_label_transfer_tbl",
                                  "doublet_identification_tbl"),
      other_arguments_to_map = c("empty_droplets_tbl",
                                 "non_batch_variation_removal_S",
                                 "alive_identification_tbl",
                                 "cell_cycle_score_tbl",
                                 "annotation_label_transfer_tbl",
                                 "doublet_identification_tbl")
    ) 
    
  )
  
}


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
  
  if(! "remove_empty_DropletUtils" %in% names(input_hpc))
    target_chunk_undefined_remove_empty_DropletUtils(input_hpc)
  
  #-----------------------#
  # Annotate cell type
  #-----------------------#
  
  if(
    !("annotate_cell_type" %in% names(input_hpc) |
    ( "remove_dead_scuttle" %in% names(input_hpc) & !is.null(input_hpc$remove_dead_scuttle$group_by))
  ))
      target_chunk_undefined_annotate_cell_type(input_hpc)
  
  #-----------------------#
  # Remove dead
  #-----------------------#
  
  if(! "remove_dead_scuttle" %in% names(input_hpc))
      target_chunk_undefined_remove_dead_scuttle(input_hpc)
  
  
  #-----------------------#
  # score cell cycle
  #-----------------------#
  if(! "score_cell_cycle_seurat" %in% names(input_hpc))
      target_chunk_undefined_score_cell_cycle_seurat(input_hpc)
  
  #-----------------------#
  # Doublets
  #-----------------------#
  
  if(! "remove_doublets_scDblFinder" %in% names(input_hpc))
     target_chunk_undefined_remove_doublets_scDblFinder(input_hpc)
  
  #-----------------------#
  # SCT
  #-----------------------#
  
  if(! "normalise_abundance_seurat_SCT" %in% names(input_hpc))
    target_chunk_undefined_normalise_abundance_seurat_SCT(input_hpc)
  
  #-----------------------#
  # Create single cell output
  #-----------------------#
  
  tar_tier_append(
    quote(preprocessing_output_factory),
    input_hpc$initialisation$tier |> get_positions() ,
    glue("{input_hpc$initialisation$store}.R")
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
    "pseudobulk_group_by.rds",
    "temp_tiers.rds",
    "temp_gene_nomenclature.rds"
  ) |> 
    remove_files_safely()
  
  return(
    tar_meta(store = glue("{input_hpc$initialisation$store}")) |> 
      filter(name |> str_detect("preprocessing_output_S_?.*$")) 
  )
}

#' @importFrom methods show
#' @export
print.HPCell <- function(x, ...){
  
  x |>
    evaluate_hpc() |> 
    print()
}
