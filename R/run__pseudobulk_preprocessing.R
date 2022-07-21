
# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_path_preprocessing_output = args[2:(length(args)-1)]
output_path = args[[length(args)]]

renv::load(project = code_directory)

library(tidyverse)
library(Seurat)
library(glue)
library(tidyseurat)
library(tidySingleCellExperiment)
library(tidysc)
library(tidybulk)

# Parallelise internally
library(furrr)
options(future.globals.maxSize = 49 * 1024 ^ 3) # 50Gb
plan(strategy = "multisession", workers = 10)


# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

counts =
  input_path_preprocessing_output |>

  # Aggregate
  future_map_dfr(~
        readRDS(.x) |>
          aggregate_cells(c(sample, predicted.celltype.l2), slot = "counts", assays=c("RNA", "ADT"))
    ) |>

  # Reshape to make RNA and ADT both features
  pivot_longer(
    cols = c(abundance_ADT, abundance_RNA),
    names_to = "data_source",
    values_to = "count"
  ) |>
  filter(!count |> is.na()) |>

  # Some manipulation to get unique feature because RNa and ADT both can have sma name genes
  rename(symbol = transcript) |>
  mutate(data_source = str_remove(data_source, "abundance_")) |>
  unite( ".feature", c(symbol, data_source), remove = FALSE)

# cell_types = counts |> distinct(predicted.celltype.l2) |> pull(predicted.celltype.l2) # |> rlang::syms()
#
# # Create SummarizedExperiment
# counts =
#   counts |>
#
#   # Spread and convert
#   select(-.aggregated_cells) |>
#   select(-predicted.celltype.l1) |>
#   pivot_wider(names_from = predicted.celltype.l2, values_from = count)
#
# cell_types_position = colnames(counts) %in% cell_types |> which()

counts |>
  as_SummarizedExperiment(
    .sample = c( sample,predicted.celltype.l2),
    .transcript = .feature,
    .abundance = count
  ) |>

  # Save
  saveRDS(output_path)

