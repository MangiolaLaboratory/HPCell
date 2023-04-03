# Rscript analysis/split_single_cell_object_into_files.R data/alive_identification/counts_alive.rds data/alive_identification _alive

# This script slipts the dataset and creates the list of files in a specific directory

set.seed(42)

# Read arguments
# args = commandArgs(trailingOnly=TRUE)
# code_directory = args[[1]]
# input_path_demultiplexed = args[[2]]
# input_path_empty_droplets = args[[3]]
# input_path_annotation_label_transfer = args[[4]]
# output_path = args[[5]]
# 
# renv::load(project = code_directory)

library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(dplyr); library(tidyr); library(ggplot2)
library(purrr)
library(Seurat)
library(tidyseurat)
library(glue)
library(scater)
library(stringr)


alive_identification<-function(code_directory,input_path_demultiplexed,input_path_empty_droplets,input_path_annotation_label_transfer,output_path){
  # Create dir
  output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)
  
  input_file =
    readRDS(input_path_demultiplexed) |>
    left_join(readRDS(input_path_empty_droplets), by=".cell") |>
    tidyseurat::filter(!empty_droplet)
  
  location <- mapIds(
    EnsDb.Hsapiens.v86,
    keys=rownames(input_file),
    column="SEQNAME",
    keytype="SYMBOL"
  )
  
  which_mito = rownames(input_file) |> str_which("^MT-")
  
  mitochondrion =
    input_file |>
    GetAssayData( slot = "counts", assay="RNA") |>
    
    # Join mitochondrion statistics
    scuttle::perCellQCMetrics(subsets=list(Mito=which_mito)) |>
    as_tibble(rownames = ".cell") |>
    dplyr::select(-sum, -detected) |>
    
    # Join cell types
    left_join(readRDS(input_path_annotation_label_transfer), by = ".cell") |>
    
    # Label cells
    nest(data = -blueprint_first.labels.fine) |>
    mutate(data = map(
      data,
      ~ .x |>
        mutate(high_mitochondrion = isOutlier(subsets_Mito_percent, type="higher")) |>
        
        # For compatibility
        mutate(high_mitochondrion = as.logical(high_mitochondrion))
    )) |>
    unnest(data)
  
  
  ribosome =
    
    input_file |>
    
    tidyseurat::select(.cell) |>
    
    # Join mitochondrion statistics
    mutate(mito_RPS = PercentageFeatureSet(input_file,  pattern = "^RPS|^RPL", assay = "RNA")[,1]  ) |>
    
    # Join cell types
    left_join(readRDS(input_path_annotation_label_transfer), by = ".cell") |>
    
    # Label cells
    nest(data = -blueprint_first.labels.fine) |>
    mutate(data = map(
      data,
      ~ .x |>
        # Label cells
        mutate(high_RPS = isOutlier(mito_RPS, type="higher")) |>
        
        # For compatibility
        mutate(high_RPS = as.logical(high_RPS)) |>
        
        as_tibble() |>
        dplyr::select(.cell, mito_RPS, high_RPS)
    )) |>
    unnest(data)
  
  
  
  
  # Save
  mitochondrion |>
    left_join(ribosome, by=".cell") |>
    mutate(alive = !high_mitochondrion & !high_RPS ) |>
    saveRDS(output_path)
  
  return(output_path)
}


