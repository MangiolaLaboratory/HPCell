set.seed(42)



# renv::load(project = code_directory)

library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(dplyr); library(tidyr); library(ggplot2)
library(purrr)
library(Seurat)
library(tidyseurat)
library(glue)
library(scater)

non_batch_variation_removal<-function(code_directory,input_path_demultiplexed, input_path_empty_droplets,input_path_alive,input_cell_cycle_scoring,input_path_marged_variable_genes,output_path){
  # Create dir
  output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)
  
  counts =
    readRDS(input_path_demultiplexed) |>
    left_join(readRDS(input_path_empty_droplets), by = ".cell") |>
    tidyseurat::filter(!empty_droplet) |>
    
    left_join(
      readRDS(input_path_alive) |>
        select(.cell, subsets_Ribo_percent, subsets_Mito_percent),
      by=".cell"
    ) |>
    
    left_join(
      readRDS(input_cell_cycle_scoring) |>
        select(.cell, G2M.Score),
      by=".cell"
    )
  
  # tidyseurat::filter(!high_mitochondrion | !high_ribosome)
  
  variable_features = readRDS(input_path_marged_variable_genes)
  
  # Set variable features
  VariableFeatures(counts) = variable_features
  
  counts |>
    
    # Normalise RNA
    SCTransform(
      assay="RNA",
      return.only.var.genes=FALSE,
      #residual.features = variable_features,
      vars.to.regress = c("subsets_Mito_percent", "subsets_Ribo_percent", "G2M.Score"),
      vst.flavor = "v2",
      
      scale_factor=2186
    ) |>
    
    # Normalise antibodies
    when(
      "ADT" %in% names(.@assays) ~ NormalizeData(., normalization.method = 'CLR', margin = 2, assay="ADT") ,
      ~ (.)
    ) |>
    
    # Drop alive columns
    select(-subsets_Ribo_percent, -subsets_Mito_percent, -G2M.Score) |>
    
    # Save
    saveRDS(output_path)
}