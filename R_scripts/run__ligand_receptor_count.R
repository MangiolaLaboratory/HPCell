
# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_path_preprocessing = args[[2]]
output_path = args[[3]]

# renv::load(project = code_directory)

library(dplyr); library(tidyr); library(ggplot2); library(purrr)
library(Seurat)
library(glue)
library(CellChat)
library(tidyseurat)
library(tidySingleCellExperiment)
library(stringr)

#library(furrr)
#plan(multisession, workers=3)

# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)


counts_cellchat =
  input_path_preprocessing |>
  readRDS() |>

  # Count
  ligand_receptor_count(counts, predicted.celltype.l2, assay = "SCT") |>

  # Save
  saveRDS(output_path)

