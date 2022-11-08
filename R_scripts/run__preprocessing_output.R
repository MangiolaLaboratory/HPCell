
# Read arguments
args = commandArgs(trailingOnly = TRUE)
code_directory = args[[1]]
input_path_non_batch_variation_removal = args[[2]]
input_path_alive = args[[3]]
input_path_annotation_label_transfer = args[[4]]
input_path_doublet_identification = args[[5]]
output_path = args[[6]]

#renv::load(project = code_directory)

library(dplyr); library(tidyr); library(ggplot2)
library(Seurat)
library(glue)
library(scDblFinder)
library(tidyseurat)
library(tidySingleCellExperiment)

# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

readRDS(input_path_non_batch_variation_removal) |>

  # Filter dead cells
  left_join(readRDS(input_path_alive) |> select(.cell, alive), by = ".cell") |>
  filter(alive) |>

  # Filter doublets
  left_join(readRDS(input_path_doublet_identification) |> select(.cell, scDblFinder.class), by = ".cell") |>
  filter(scDblFinder.class=="singlet") |>

  # Filter Red blood cells and platelets
  left_join(readRDS(input_path_annotation_label_transfer), by = ".cell") |>
  filter(!predicted.celltype.l2 %in% c("Eryth", "Platelet")) |>

  # Save
  saveRDS(output_path)

