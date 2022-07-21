
# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_path_non_batch_variation_removal = args[[2]]
input_path_alive = args[[3]]
input_path_annotation_label_transfer = args[[4]]
input_path_doublet_identification = args[[5]]
output_path = args[[6]]

renv::load(project = code_directory)

library(tidyverse)
library(Seurat)
library(glue)
library(scDblFinder)
library(tidyseurat)
library(tidySingleCellExperiment)

# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

preprocessing_output =
  readRDS(input_path_non_batch_variation_removal) |>

  # Filter dead
  left_join(
    readRDS(input_path_alive) |>

      # Otherwise the joining down the pipeline fails because the variables are of the strange class otlr.flt
      mutate(
        high_mitochondrion = as.logical(high_mitochondrion),
        high_RPS = as.logical(high_RPS)
      ),
    by = ".cell"
  ) |>
  filter(!high_mitochondrion & !high_RPS) |>

  # Annotate
  left_join(readRDS(input_path_annotation_label_transfer), by = ".cell") |>

  # Filter doublets
  filter(.cell %in% (
    readRDS(input_path_doublet_identification) |>
      filter( scDblFinder.class == "singlet") |>
      pull(.cell)
  ))

preprocessing_output |>

  # Save
  saveRDS(output_path)

