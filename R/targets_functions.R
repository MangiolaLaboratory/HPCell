
#' @importFrom targets tar_script
#' @importFrom targets tar_option_set
#' @importFrom dplyr pull
#' @importFrom dplyr count
#' @importFrom dplyr rename
#' 
#' @import targets
#' 
#' @export
#' 
hpcell_map_test_differential_abundance = function(data_df, formula, .data_column, store =  tempfile(tmpdir = "."), append = FALSE){
  
  .data_column = enquo(.data_column)
  
  # Check that names are different
  if(data_df |> count(name) |> pull(n) |> max() > 1)
    stop("HPCell says: the column name must contain unique identifiers")
  
  data_df |> rename(data = !!.data_column) |>  saveRDS("temp_data.rds")
  formula |>  saveRDS("temp_formula.rds")
  

  # library(future)
  # library("future.batchtools")
  # 
  # small_slurm = 
  #   tar_resources(
  #     future = tar_resources_future(
  #       plan = tweak(
  #         batchtools_slurm,
  #         template = glue("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/third_party_sofware/slurm_batchtools.tmpl"),
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
  #         template = glue("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/third_party_sofware/slurm_batchtools.tmpl"),
  #         resources = list(
  #           ncpus = 19,
  #           memory = 6000,
  #           walltime = 172800
  #         )
  #       )
  #     )
  #   )
  
  
  # Header
  if(!append)
    tar_script_append({
    
    
    #-----------------------#
    # Input
    #-----------------------#
    library(targets)
    library(tarchetypes)
    
    #-----------------------#
    # Future SLURM
    #-----------------------#
    
    library(future)
    library(future.callr)
    plan(callr)
    # library("future.batchtools")
    # slurm <- 
    # 	`batchtools_slurm` |>
    # 	future::tweak( template = glue("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/third_party_sofware/slurm_batchtools.tmpl"),
    # 								 resources=list(
    # 								 	ncpus = 20,
    # 								 	memory = 6000,
    # 								 	walltime = 172800
    # 								 )
    # 	)
    # plan(slurm)
    
    #-----------------------#
    # Packages
    #-----------------------#
    tar_option_set( 
      packages = c(
        "stringr", "tibble", "tidySingleCellExperiment", "dplyr", "Matrix",
        "Seurat", "tidyseurat", "glue", "purrr", "tidybulk", "tidySummarizedExperiment", "edgeR",
        "digest", "jascap"
      ), 
      storage = "worker", 
      retrieval = "worker", 
      error = "continue", 		
      format = "qs"
      #,
      #debug = "estimates_b84ea256", # Set the target you want to debug.
      #cue = tar_cue(mode = "never") # Force skip non-debugging outdated targets.
    )

    list_of_tar_de = 
      list(
        tar_target(file, "temp_data.rds", format = "file"),
        tar_target(formula, readRDS("temp_formula.rds"))
      )
    
  }, glue("{store}.R"))
  
 
  tar_script_append({
    
    #-----------------------#
    # Pipeline
    #-----------------------#
   list_of_tar_de = list_of_tar_de |> c(list(
      
      tarchetypes::tar_group_by(pseudobulk_df_tissue, readRDS(file), name),
      
      # Dispersion
      tar_target(
        pseudobulk_df_tissue_dispersion, 
        pseudobulk_df_tissue |> map_add_dispersion_to_se(data), 
        pattern = map(pseudobulk_df_tissue),
        iteration = "group"
        #, 
        #resources = small_slurm
      ),
      
      # Split in gene chunks
      tar_target(
        pseudobulk_df_tissue_split_by_gene, 
        pseudobulk_df_tissue_dispersion |> map_split_se_by_gene(data), 
        pattern = map(pseudobulk_df_tissue_dispersion),
        iteration = "group"
        #, 
        #resources = small_slurm
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
            max_rows_for_matrix_multiplication = 10000, 
            cores = 18
          ) , 
        pattern = map(pseudobulk_df_tissue_split_by_gene_grouped),
        iteration = "group"
        #, 
        #resources = big_slurm
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
        #, 
        #resources = big_slurm
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
  if(!append)
  tar_make_future(
    script = glue("{store}.R"),
    store = store, 
    workers = 10
    #, 
    #workers = 200, 
    #garbage_collection = TRUE
  )
  
  if(!append)
  file.remove("temp_data.rds")
  file.remove("temp_formula.rds")
  
  if(!append)
  tar_read(estimates, store = store)
}

#' @importFrom tibble tibble
#' 
#' @export
#' 
hpcell_test_differential_abundance = function(.data, formula, store =  tempfile(tmpdir = "."), append = FALSE){
  
  # Create input dataframe
  tibble(name = "my_data", data = list(!!.data )) |> 
    
    # Call map function 
    hpcell_map_test_differential_abundance(formula, data, store = store, append = append)
  
}


