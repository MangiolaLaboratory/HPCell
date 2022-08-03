renv::activate()

library(tidyverse)
library(Seurat)
library(tidyseurat)
library(glue)

set.seed(42)


# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_path_demultiplexed = args[[2]]
input_path_empty_droplets = args[[3]]
cell_type_column_for_subsetting = args[[4]]
output_path = args[[5]]

# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

# Read data and filter
readRDS(input_path_demultiplexed) |>
  left_join(readRDS(input_path_empty_droplets), by = ".cell") |>
  tidyseurat::filter(!empty_droplet)


input_file = readRDS(input_path)

input_file %>%

  # Sample up to a plateau to avoid extreme cell_type bias
  nest(data = -!!as.symbol(cell_type_column_for_subsetting)) %>%
  mutate(n = map_int(data, ~ ncol(.x))) %>%
  mutate(upper_quantile = quantile(n, 0.75) %>% as.integer()) %>%
  mutate(data = map2(
    data, upper_quantile,
    ~ sample_n(.x, min(ncol(.x), .y), replace = FALSE)
  )) %>%
  filter(n>1) %>%
  unnest(data) %>%

  # Normalise
  NormalizeData(assay="RNA") %>%

  # Get variable features
  FindVariableFeatures(assay="RNA", nfeatures = 1000) %>%
  VariableFeatures(assay = "RNA") %>%
  saveRDS(output_path)
