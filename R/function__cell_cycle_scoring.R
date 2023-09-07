# Rscript analysis/split_single_cell_object_into_files.R data/alive_identification/counts_alive.rds data/alive_identification _alive

# library(future)
# options(future.globals.maxSize = 30 * 1024 ^ 3) # 30Gb
# plan(strategy = "multisession", workers = 10)


# renv::load(project = code_directory)

library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(dplyr); library(tidyr); library(ggplot2)
library(purrr)
library(Seurat)
library(tidyseurat)
library(glue)
library(scater)

cell_cycle_scoring<-function(code_directory, input_path_demultiplexed, input_path_empty_droplets,output_path){

  # Create dir
  output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)
  
  counts =
    readRDS(input_path_demultiplexed) |>
    left_join(readRDS(input_path_empty_droplets), by = ".cell") |>
    tidyseurat::filter(!empty_droplet) |>
    
    # Normalise needed
    NormalizeData() |>
    
    # Scoring
    CellCycleScoring(
      s.features = cc.genes$s.genes,
      g2m.features = cc.genes$g2m.genes,
      set.ident = FALSE
    ) |>
    
    as_tibble() |>
    select(.cell,  S.Score, G2M.Score, Phase) |>
    
    # Save
    saveRDS(output_path)

}