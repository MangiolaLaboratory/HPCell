set.seed(42)

# Read arguments
args = commandArgs(trailingOnly = TRUE)
code_directory = args[[1]]
input_path_demultiplexed = args[[2]]
input_path_empty_droplets = args[[3]]
input_path_marged_variable_genes = args[[4]]
output_path = args[[5]]

renv::load(project = code_directory)

library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(dplyr); library(tidyr); library(ggplot2)
library(purrr)
library(Seurat)
library(tidyseurat)
library(glue)
library(scater)

code_directory = args[[1]]
input_path_demultiplexed = args[[2]]
input_path_empty_droplets = args[[3]]
input_path_marged_variable_genes = args[[4]]
output_path = args[[5]]
non_batch_variation_removal<-function(code_directory,tissue, input_paths,cell_type_column_for_subsetting,output_path){
  # Create dir
  output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)
  
  counts =
    readRDS(input_path_demultiplexed) |>
    left_join(readRDS(input_path_empty_droplets), by = ".cell") |>
    tidyseurat::filter(!empty_droplet)
  
  # left_join(readRDS(input_path_alive), by=".cell") |>
  # tidyseurat::filter(!high_mitochondrion | !high_RPS)
  
  variable_features = readRDS(input_path_marged_variable_genes)
  
  # Set variable features
  VariableFeatures(counts) = variable_features
  
  counts |>
    
    # Normalise RNA
    SCTransform(assay="RNA", return.only.var.genes=FALSE, residual.features = variable_features) |>
    
    # Normalise antibodies
    when(
      "ADT" %in% names(.@assays) ~ NormalizeData(., normalization.method = 'CLR', margin = 2, assay="ADT") ,
      ~ (.)
    ) |>
    
    # Save
    saveRDS(output_path)
}