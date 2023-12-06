
#' Main Function for HPCell Map Test Differential Abundance
#'
#' @description
#' This function prepares and runs a differential abundance test pipeline using the 'targets' package. It sets up necessary files, appends scripts, and executes the pipeline.
#'
#' @param data_df Data frame to be processed.
#' @param formula Formula for the differential abundance test.
#' @param .data_column Column in the data frame containing the data.
#' @param store File path for temporary storage.
#' @param computing_resources Computing resources configuration.
#' @param cpus_per_task Number of CPUs allocated per task.
#' @param debug_job_id Optional job ID for debugging.
#' @param append Flag to append to existing script.
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
#' 
#' @export
hpcell_map_test_differential_abundance = function(
    data_df,
    formula, 
    .data_column, 
    .group_name_columns,
    .abundance = NULL,
    store =  tempfile(tmpdir = "."), 
    computing_resources = crew_controller_local(workers = 1) , 
    cpus_per_task = 1,
    debug_job_id = NULL,
    append = FALSE
  ){
  
  .data_column = enquo(.data_column)
  .group_name_columns = enquo(.group_name_columns)
  .abundance = enquo(.abundance)
  
  if(quo_is_symbolic(.abundance)) .abundance = quo_names(.abundance)
  else .abundance =  
    data_df |> 
    pull(data) |> 
    extract2(1) |> 
    assays() |> 
    names() |> 
    extract2(1)
  
  # Check that names are different
  if(data_df |> count(!!.group_name_columns) |> pull(n) |> max() > 1)
    stop("HPCell says: the column name must contain unique identifiers")
  
  data_df |> rename(data = !!.data_column, name = !!.group_name_columns) |>  saveRDS("temp_data.rds")
  formula |>  saveRDS("temp_formula.rds")
  computing_resources |> saveRDS("temp_computing_resources.rds")
  debug_job_id |> saveRDS("temp_debug_job_id.rds")
  .abundance |> saveRDS("temp_abundance_column_name.rds")
  



  
  # Header
  if(!append)
    tar_script_append({
    
    
    #-----------------------#
    # Input
    #-----------------------#
    library(targets)
    library(tarchetypes)
    
    computing_resources = readRDS("temp_computing_resources.rds")
    debug_job_id = readRDS("temp_debug_job_id.rds")
    
    #-----------------------#
    # Packages
    #-----------------------#
    tar_option_set( 
      packages = c(
        "stringr", "tibble", "tidySingleCellExperiment", "dplyr", "Matrix",
        "Seurat", "tidyseurat", "glue", "purrr", "tidybulk", "tidySummarizedExperiment", "edgeR",
        "digest", "HPCell"
      ), 
      storage = "worker", 
      retrieval = "worker", 
      # error = "continue", 		
      format = "qs",
      controller = computing_resources,
      debug = debug_job_id # Set the target you want to debug.
      #cue = tar_cue(mode = "never") # Force skip non-debugging outdated targets.
    )

    list_of_tar_de = 
      list(
        tar_target(file, "temp_data.rds", format = "file"),
        tar_target(formula, readRDS("temp_formula.rds")),
        tar_target(abundance, readRDS("temp_abundance_column_name.rds"))
      )
    
  }, glue("{store}.R"))
  
 
  tar_script_append({
    
    #-----------------------#
    # Pipeline
    #-----------------------#
   list_of_tar_de = list_of_tar_de |> c(list(
      
      tarchetypes::tar_group_by(pseudobulk_df_tissue, readRDS(file), name),
      
      tar_target( computing_resources, readRDS("temp_computing_resources.rds") ),
      
      # Dispersion
      tar_target(
        pseudobulk_df_tissue_dispersion, 
        pseudobulk_df_tissue |> map_add_dispersion_to_se(data, abundance), 
        pattern = map(pseudobulk_df_tissue),
        iteration = "group"
        # , 
        # resources = computing_resources
      ),
      
      # Split in gene chunks
      tar_target(
        pseudobulk_df_tissue_split_by_gene, 
        pseudobulk_df_tissue_dispersion |> map_split_se_by_gene(
          data, 
          computing_resources$client$workers
        ), 

        pattern = map(pseudobulk_df_tissue_dispersion),
        iteration = "group"
        # , 
        # resources = computing_resources
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
          map_test_differential_abundance(
            data,
            formula, 
            .abundance = abundance,
            max_rows_for_matrix_multiplication = 10000, 
            cores = 1
          ) , 
        pattern = map(pseudobulk_df_tissue_split_by_gene_grouped),
        iteration = "group"
        # , 
        # resources = computing_resources
      ),
      
      # Regroup
      tarchetypes::tar_group_by(
        estimates_regrouped, 
        estimates_chunk, 
        name
      ),
      tar_target(
        estimates, 
        estimates_regrouped |> unnest_summarized_experiment(data) |> nest(se = -name) , 
        pattern = map(estimates_regrouped),
        iteration = "group"
        # , 
        # resources = computing_resources
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
  
  if(!append)
  file.remove("temp_data.rds")
  file.remove("temp_formula.rds")
  file.remove("temp_computing_resources.rds")
  
  
  if(!append)
  tar_read(estimates, store = store)
}

#' Wrapper Function for HPCell Test Differential Abundance
#'
#' @description
#' A wrapper function that formats data into a tibble and calls `hpcell_map_test_differential_abundance` for differential abundance testing.
#'
#' @param .data Data frame or similar object for analysis.
#' @param formula Formula for the differential abundance test.
#' @param store File path for temporary storage.
#' @param computing_resources Computing resources configuration.
#' @param cpus_per_task Number of CPUs allocated per task.
#' @param debug_job_id Optional job ID for debugging.
#' @param append Flag to append to existing script.
#'
#' @return A `targets` pipeline output, typically a nested tibble with differential abundance estimates.
#'
#' @importFrom tibble tibble
#' @export
#' 
hpcell_test_differential_abundance = function(
    .data, 
    formula,
    store =  tempfile(tmpdir = "."),  
    computing_resources = crew_controller_local(workers = 1) , 
    debug_job_id = NULL, 
    append = FALSE
  ){
  
  # Create input dataframe
  tibble(name = "my_data", data = list(!!.data )) |> 
    
    # Call map function 
    hpcell_map_test_differential_abundance(
      formula,
      data,    
      .group_name_columns = name,
      .abundance = NULL, 
      store = store, 
      computing_resources = computing_resources,
      debug_job_id = debug_job_id, 
      append = append
    )

  
}


