# Rscript analysis/split_single_cell_object_into_files.R data/alive_identification/counts_alive.rds data/alive_identification _alive

# library(future)
# options(future.globals.maxSize = 30 * 1024 ^ 3) # 30Gb
# plan(strategy = "multisession", workers = 10)

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
library(tidyverse)
library(purrr)
library(Seurat)
library(tidyseurat)
library(glue)
library(scater)

# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

readRDS(input_path_demultiplexed) |>
  left_join(readRDS(input_path_empty_droplets), by = ".cell") |>
  tidyseurat::filter(!empty_droplet) |>

  # left_join(readRDS(input_path_alive), by=".cell") |>
  # tidyseurat::filter(!high_mitochondrion | !high_RPS) |>

  # Normalise RNA
  SCTransform(assay="RNA", residual.features = readRDS(input_path_marged_variable_genes), return.only.var.genes=FALSE) |>

  # Normalise antibodies
  when(
    "ADT" %in% names(.@assays) ~ NormalizeData(., normalization.method = 'CLR', margin = 2, assay="ADT") ,
    ~ (.)
  ) |>

  # Save
  saveRDS(output_path)


