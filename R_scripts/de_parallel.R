library(tidyverse)
library(targets)
library(glue)


# Get input

my_data = 
  tidySummarizedExperiment::se |> 
  tidybulk::keep_abundant()
tibble(name = "my_data", se = list(my_data)) |> saveRDS(glue("temp_data.rds"))
my_store = tempfile(tmpdir = ".")

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

tar_script({
  
  
  #-----------------------#
  # Input
  #-----------------------#
  library(tidyverse)
  library(targets)
  library(tarchetypes)
  library(tidyseurat)
  library(glue)
  library(qs)
  
  

  
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
      "digest"
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
  # Functions
  #-----------------------#
  se_add_dispersion = function(se_df){
    
    se_df |> 
      mutate(se = map(
        se,
        ~ {
          counts = .x |> assay("counts")
          
          .x |> 
            left_join(
              
              # Dispersion data frame
              estimateDisp(counts)$tagwise.dispersion |> 
                setNames(rownames(counts)) |>
                enframe(name = ".feature", value = "dispersion")
            )
        }
      ))
    
  }
  
  split_by_gene = function(se_df){
    se_df |> 
      mutate(se = map(
        se,
        ~ {
          chunks =
            tibble(.feature = rownames(.x)) |>
            mutate(chunk___ = sample(1:10, n(), replace = TRUE))
          
          .x |> 
           left_join(chunks) |>
          nest(se_chunk = -chunk___)
        }
      )) |>
      unnest(se) |> 
      select(-chunk___) |>
      mutate(se_md5 = map_chr(se_chunk, digest))
  }
  
  analyse = function(
    se, max_rows_for_matrix_multiplication = NULL, 
    cores = 1
  ){
    
    se |> mutate(se_chunk = map(
      se_chunk,
      ~ .x |>
        
        # Test
        test_differential_abundance(
          ~ dex + (1 | cell),
          method = "glmmSeq_lme4",
          cores = cores, 
          max_rows_for_matrix_multiplication = max_rows_for_matrix_multiplication,
          .dispersion = dispersion
        )
      
    ))
    
 
    
     
    
    
  }
  
  #-----------------------#
  # Pipeline
  #-----------------------#
  list(
    
    
    tar_target(file, "temp_data.rds", format = "file"),
    tarchetypes::tar_group_by(pseudobulk_df_tissue, readRDS(file), name),
  
    # Dispersion
    tar_target(
      pseudobulk_df_tissue_dispersion, 
      pseudobulk_df_tissue |> se_add_dispersion(), 
      pattern = map(pseudobulk_df_tissue)
      #, 
      #resources = small_slurm
    ),
    
    # Split in gene chunks
    tar_target(
      pseudobulk_df_tissue_split_by_gene, 
      pseudobulk_df_tissue_dispersion |> split_by_gene(), 
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
      estimates, 
      pseudobulk_df_tissue_split_by_gene_grouped |> analyse(max_rows_for_matrix_multiplication = 10000, cores = 18) , 
      pattern = map(pseudobulk_df_tissue_split_by_gene_grouped)
      #, 
      #resources = big_slurm
    )
  )
  
  
}, ask = FALSE, script = glue("{my_store}.R"))


tar_make_future(
  script = glue("{my_store}.R"),
  store = my_store, 
  workers = 10
  #, 
  #workers = 200, 
  #garbage_collection = TRUE
)
