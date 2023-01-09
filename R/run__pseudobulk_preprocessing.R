
# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_path_preprocessing_output = args[2:(length(args)-2)]
output_path_sample_cell_type = args[[length(args)-1]]
output_path_sample = args[[length(args)]]

renv::load(project = code_directory)

library(dplyr); library(tidyr); library(ggplot2)
library(Seurat)
library(glue)
library(tidyseurat)
library(tidySingleCellExperiment)
library(tidysc)
library(tidybulk)
library(rlang)
library(stringr)
library(purrr)

# # Parallelise internally
# library(furrr)
# options(future.globals.maxSize = 49 * 1024 ^ 3) # 50Gb
# plan(strategy = "multisession", workers = 10)

# Do we have RNA and also ADT
assays = readRDS(input_path_preprocessing_output[1])@assays |> names() |> intersect(c("RNA", "ADT"))

# Create dir
output_path_sample |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)


input_path_preprocessing_output |>

  # Aggregate
  map_dfr(~ {
    library(rlang)
    readRDS(.x) |>
      aggregate_cells(c(sample, predicted.celltype.l2), slot = "counts", assays=assays)
  }) |>

  # Reshape to make RNA and ADT both features
  pivot_longer(
    cols = assays,
    names_to = "data_source",
    values_to = "count"
  ) |>
  filter(!count |> is.na()) |>

  # Some manipulation to get unique feature because RNa and ADT both can have sma name genes
  rename(symbol = feature) |>
  mutate(data_source = str_remove(data_source, "abundance_")) |>
  unite( ".feature", c(symbol, data_source), remove = FALSE) |>

  # Covert
  as_SummarizedExperiment(
    .sample = c( sample,predicted.celltype.l2),
    .transcript = .feature,
    .abundance = count
  ) |>

  # Save
  saveRDS(output_path_sample_cell_type)

gc()

# ONLY SAMPLE
input_path_preprocessing_output |>

  # Aggregate
  map_dfr(~
                   readRDS(.x) |>
                   aggregate_cells(c(sample), slot = "counts", assays=assays)
  ) |>

  # Reshape to make RNA and ADT both features
  pivot_longer(
    cols = assays,
    names_to = "data_source",
    values_to = "count"
  ) |>
  filter(!count |> is.na()) |>

  # Some manipulation to get unique feature because RNa and ADT both can have sma name genes
  rename(symbol = feature) |>
  mutate(data_source = str_remove(data_source, "abundance_")) |>
  unite( ".feature", c(symbol, data_source), remove = FALSE) |>

  # Covert
  as_SummarizedExperiment(
    .sample = c( sample),
    .transcript = .feature,
    .abundance = count
  ) |>

  # Save
  saveRDS(output_path_sample)

