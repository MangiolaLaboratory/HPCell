
#' Main Function for HPCell Map Test Differential Abundance
#'
#' @description
#' This function prepares and runs a differential abundance test pipeline using the 'targets' package. It sets up necessary files, appends scripts, and executes the pipeline.
#'
#' @param formula_list List of formula for the differential abundance test.
#' @param store File path for temporary storage.
#' @param computing_resources Computing resources configuration.
#' @param cpus_per_task Number of CPUs allocated per task.
#' @param debug_job_id Optional job ID for debugging.
#' @param append Flag to append to existing script.
#' @param data_list list of dataframes to be processed
#' @param .abundance (optional) A symbol or string indicating the column name in the `SingleCellExperiment` object to be used for abundance measures. If not explicitly provided, the function attempts to automatically detect an appropriate column by examining the first object in `data_list`.
#' @param ... additional arguments 
#'
#' @return A `targets` pipeline output, typically a nested tibble with differential abundance estimates.
#'
#' @importFrom targets tar_script
#' @importFrom targets tar_option_set
#' @importFrom dplyr pull
#' @importFrom dplyr count
#' @importFrom dplyr rename
#' @importFrom crew crew_controller_local
#' @importFrom magrittr extract2
#' @import targets
#' @importFrom rlang quo_is_symbolic
#' @importFrom SummarizedExperiment assays
#' @importFrom tibble rowid_to_column
#' @importFrom callr r
#' @importFrom tibble rowid_to_column
#' 
#' @export
map2_test_differential_abundance_hpc = function(
    data_list,
    formula_list, 
    .abundance = NULL,
    store =  tempfile(tmpdir = "."), 
    computing_resources = crew_controller_local(workers = 1) , 
    cpus_per_task = 1,
    debug_job_id = NULL,
    append = FALSE,
    ...
  ){
  
  #Fix GChecks 
  abundance = NULL 
  file_data = NULL 
  file_formula = NULL
  .abundance = NULL
  number_of_workers = NULL
  number_of_datasets = NULL
  pseudobulk_df_tissue = NULL
  name = NULL 
  pseudobulk_df_tissue_dispersion = NULL 
  pseudobulk_df_tissue_split_by_gene = NULL
  pseudobulk_df_tissue_split_by_gene_grouped = NULL
  se_md5 = NULL
  estimates_chunk = NULL
  my_group = NULL
  assay_name = NULL 
  
  .abundance = enquo(.abundance)

  if(quo_is_symbolic(.abundance)) .abundance = quo_names(.abundance)
  else .abundance =  
    data_list[[1]] |> 
    assays() |> 
    names() |> 
    extract2(1)
  
  data_list |> saveRDS("temp_data.rds")
  
  # convert to character because formula captures some 
  # of the local environment and creates very big files
  formula_list |> map(as.character) |>  saveRDS("temp_formula.rds")
  computing_resources |> saveRDS("temp_computing_resources.rds")
  debug_job_id |> saveRDS("temp_debug_job_id.rds")
  .abundance |> saveRDS("temp_abundance_column_name.rds")
  data_list |> length() |> saveRDS("temp_number_of_datasets.rds")
  

  # Header
  if(!append)
    tar_script_append({
    
    
    #-----------------------#
    # Input
    #-----------------------#
    # library(targets)
    # library(tarchetypes)
    
    computing_resources = readRDS("temp_computing_resources.rds")
    debug_job_id = readRDS("temp_debug_job_id.rds")
    
    #-----------------------#
    # Packages
    #-----------------------#
    tar_option_set( 
      packages = c(
        "stringr", "tibble", "tidySingleCellExperiment", "dplyr", "tidyseurat", "glue", "purrr", "tidybulk", "tidySummarizedExperiment", "edgeR",
        "digest", "HPCell"
      ), 
      storage = "worker", 
      retrieval = "worker", 
      # error = "continue", 		
      format = "qs",
      controller = computing_resources,
      resources = tar_resources(
        qs = tar_resources_qs(preset = "fast")
      ),
      debug = debug_job_id # Set the target you want to debug.
      #cue = tar_cue(mode = "never") # Force skip non-debugging outdated targets.
    )

    list_of_tar_de = 
      list(
        tar_target(file_data, "temp_data.rds", format = "file", deployment = "main"),
        tar_target(file_formula, "temp_formula.rds", format = "file", deployment = "main"),
        tar_target(abundance, readRDS("temp_abundance_column_name.rds"), deployment = "main"),
        tar_target( number_of_workers, readRDS("temp_computing_resources.rds")$client$workers, deployment = "main" ),
        tar_target( number_of_datasets, readRDS("temp_number_of_datasets.rds"), deployment = "main" )
      )
    
  }, glue("{store}.R"))
  
 
  tar_script_append({
    
    #-----------------------#
    # Pipeline
    #-----------------------#
   list_of_tar_de = list_of_tar_de |> c(list(
      
      tarchetypes::tar_group_by(
        pseudobulk_df_tissue, 
         tibble(
          data = readRDS(file_data),
          formula = readRDS(file_formula)
        ) |> 
          rowid_to_column(var = "name"), 
        name, 
        deployment = "main"
      ),
      
      # Dispersion
      tar_target(
        pseudobulk_df_tissue_dispersion, 
        pseudobulk_df_tissue |> map_add_dispersion_to_se(data, formula, abundance), 
        pattern = map(pseudobulk_df_tissue),
        iteration = "group"
      ),
      
      # Split in gene chunks
      tar_target(
        pseudobulk_df_tissue_split_by_gene, 
        pseudobulk_df_tissue_dispersion |> map_split_se_by_number_of_genes(
          data, 
          chunk_size = 100 # / number_of_datasets
        ), 

        pattern = map(pseudobulk_df_tissue_dispersion),
        iteration = "group"
      ),
      
      # Parallelise rows
      tarchetypes::tar_group_by(
        pseudobulk_df_tissue_split_by_gene_grouped, 
        pseudobulk_df_tissue_split_by_gene, 
        name, se_md5
      ),
      
      # Analyse
      tar_target(
        estimates_chunk, 
        pseudobulk_df_tissue_split_by_gene_grouped |>
          
          # transform back to formula because I converted to character before
          mutate(formula = map(formula, as.formula)) |> 
          
          map_test_differential_abundance(
            data,
            .formula = formula, 
            .abundance = abundance,
            max_rows_for_matrix_multiplication = 10000, 
            cores = 1, action = "get"
          ) |> 
          
          # For some reason it occupies a LOT of space (29Mb) 
          # probably ecause is carrying local variable with it
          select(-formula), 
        pattern = map(pseudobulk_df_tissue_split_by_gene_grouped),
        iteration = "group"
      )
      
    ))

    
  }, glue("{store}.R"))
  
  # Return the whole pipeline
  if(!append)
    tar_script_append({
    
    #-----------------------#
    # Pipeline
    #-----------------------#
    list_of_tar_de 
    
  }, glue("{store}.R"))
  
  # Execute pipeline
  if(!append){
    
     if(is.null(debug_job_id)) callr_function = callr::r
     else callr_function = NULL
    
    tar_make(
      script = glue("{store}.R"),
      store = store, 
      callr_function = callr_function
      #, 
      #workers = 200, 
      #garbage_collection = TRUE
    )
  }
  
  if(!append){
    file.remove("temp_data.rds")
    # file.remove("temp_formula.rds")
    file.remove("temp_computing_resources.rds")
  }

  message("HPCell says: Start collecting results.")

  estimates = 
    tar_read(estimates_chunk, store = store) |> 
    nest(my_group = -name) |> 
    pull(my_group) |> 
    map(~ do.call(
      rbind, 
      pull(.x, data) 
    ))
  
  # Return
  if(!append)
    map2(
      data_list,
      estimates,
      ~ {
        rowData(.x) = rowData(.x) |> cbind(.y)
        .x
      }
    )
  
}

#' Wrapper Function for HPCell Test Differential Abundance
#'
#' @description
#' A wrapper function that formats data into a tibble and calls `map2_test_differential_abundance_hpc` for differential abundance testing.
#' @importFrom magrittr extract2

#' @param .data Data frame or similar object for analysis.
#' @param formula Formula for the differential abundance test.
#' @param store File path for temporary storage.
#' @param computing_resources Computing resources configuration.
# cpus_per_task Number of CPUs allocated per task.
#' @param debug_job_id Optional job ID for debugging.
#' @param append Flag to append to existing script.
#' @param ... additional arguments 
#'
#' @return A `targets` pipeline output, typically a nested tibble with differential abundance estimates.
#'
#' @importFrom tibble tibble
#' @export
#' 
test_differential_abundance_hpc = function(
    .data, 
    formula,
    store =  tempfile(tmpdir = "."),  
    computing_resources = crew_controller_local(workers = 1) , 
    debug_job_id = NULL, 
    append = FALSE,
    ...
  ){
  
  # Create input dataframe
  tibble(
    name = "my_data", 
    data = list(!!.data ),
    formula = list(!!formula)
  ) |> 
    
    mutate(data = map2_test_differential_abundance_hpc(
      data,
      formula ,
      store = store, 
      computing_resources = computing_resources,
      debug_job_id = debug_job_id, 
      append = append,
      ...
    )) |> 
    pull(data) |> 
    extract2(1)

  
}


