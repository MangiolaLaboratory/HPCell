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
input_path_alive = args[[4]]
input_cell_cycle_scoring = args[[5]]
input_path_marged_variable_genes = args[[6]]
output_path = args[[7]]

# renv::load(project = code_directory)

library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(dplyr); library(tidyr); library(ggplot2)
library(purrr)
library(Seurat)
library(tidyseurat)
library(glue)
library(scater)

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


