
# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
tissue = args[[2]]
input_path_non_batch_variation_removal = args[[3]]
input_path_alive = args[[4]]
input_cell_cycle_scoring = args[[5]]

input_path_annotation_label_transfer = args[[6]]
input_path_doublet_identification = args[[7]]
output_path = args[[8]]

renv::load(project = code_directory)

library(dplyr); library(tidyr); library(ggplot2)
library(Seurat)
library(glue)
library(scDblFinder)
library(tidyseurat)
library(tidySingleCellExperiment)
library(purrr)
library(magrittr)

# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

readRDS(input_path_non_batch_variation_removal) |>

  # Filter dead cells
  left_join(
    readRDS(input_path_alive) |>
      select(.cell, alive, subsets_Mito_percent, subsets_Ribo_percent, high_mitochondrion, high_ribosome),
    by = ".cell"
  ) |>
  filter(alive) |>

  # Add cell cycle
  left_join(
    readRDS(input_cell_cycle_scoring),
    by=".cell"
  ) |>

  # Filter doublets
  left_join(readRDS(input_path_doublet_identification) |> select(.cell, scDblFinder.class), by = ".cell") |>
  filter(scDblFinder.class=="singlet") |>

  # Filter Red blood cells and platelets
  left_join(readRDS(input_path_annotation_label_transfer), by = ".cell") |>

  when(
    tissue |>
      tolower() |>
      equals("pbmc") ~
      filter(., !predicted.celltype.l2 %in% c("Eryth", "Platelet")),
    ~ (.)
  ) |>


  # Save
  saveRDS(output_path)

