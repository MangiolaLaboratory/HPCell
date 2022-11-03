
# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_path_demultiplexed = args[[2]]
reference_label = args[[3]]
input_path_empty_droplets = args[[4]]
input_path_alive = args[[5]]
input_path_annotation_label_transfer = args[[6]]
output_path = args[[7]]

renv::load(project = code_directory)

library(tidyverse)
library(Seurat)
library(glue)
library(scDblFinder)
library(tidyseurat)
library(tidySingleCellExperiment)

# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

readRDS(input_path_demultiplexed) |>

  # Filtering empty
  left_join(readRDS(input_path_empty_droplets), by = ".cell") |>
  filter(!empty_droplet) |>

  # Filter dead
  left_join(readRDS(input_path_alive), by = ".cell") |>
  filter(!high_mitochondrion & !high_RPS) |>

  # Annotate
  left_join(readRDS(input_path_annotation_label_transfer), by = ".cell") |>

  # Convert
  as.SingleCellExperiment() |>

  # Double identification. If no ;abel provided calculate clusters
  scDblFinder(clusters = ifelse(reference_label=="none", TRUE, reference_label)) |>

  # Parse
  as_tibble() |>
  select(.cell, contains("scDblFinder")) |>

  # Save
  saveRDS(output_path)
