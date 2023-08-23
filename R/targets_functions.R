
#' @importFrom targets tar_script
#' @importFrom targets tar_option_set
#' @import targets
#' 
#' @export
#' 
hpcell_test_differential_abundance = function(data_df, store =  tempfile(tmpdir = ".")){
  
  data_df |> saveRDS("temp_data.rds")
  
  
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
  
  substitute({
    
    
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
        "CuratedAtlasQueryR", "stringr", "tibble", "tidySingleCellExperiment", "dplyr", "Matrix",
        "Seurat", "tidyseurat", "glue", "qs",  "purrr", "tidybulk", "tidySummarizedExperiment", "edgeR",
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
    
    
    #-----------------------#
    # Pipeline
    #-----------------------#
    list(
      
      
      tar_target(file, "temp_data.rds", format = "file"),
      tarchetypes::tar_group_by(pseudobulk_df_tissue, readRDS(file), name),
      
      # Dispersion
      tar_target(
        pseudobulk_df_tissue_dispersion, 
        pseudobulk_df_tissue |> map_add_dispersion_to_se(), 
        pattern = map(pseudobulk_df_tissue)
        #, 
        #resources = small_slurm
      ),
      
      # Split in gene chunks
      tar_target(
        pseudobulk_df_tissue_split_by_gene, 
        pseudobulk_df_tissue_dispersion |> map_split_se_by_gene(), 
        pattern = map(pseudobulk_df_tissue_dispersion)
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
        pseudobulk_df_tissue_split_by_gene_grouped |> map_test_differential_abundance(max_rows_for_matrix_multiplication = 10000, cores = 18) , 
        pattern = map(pseudobulk_df_tissue_split_by_gene_grouped)
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
        estimates_regrouped |> unnest_summarized_experiment(se_chunk) |> nest(se = -name) , 
        pattern = map(estimates_regrouped)
        #, 
        #resources = big_slurm
      )
      
    )
    
    
  }) |> deparse() |> head(-1) |> tail(-1) |>  writeLines(glue("{store}.R"))
  
  
  tar_make_future(
    script = glue("{store}.R"),
    store = store, 
    workers = 10
    #, 
    #workers = 200, 
    #garbage_collection = TRUE
  )
  
  file.remove("temp_data.rds")
  
  tar_read(estimates, store = store)
}
