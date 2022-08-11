set.seed(42)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_paths = args[2:(length(args)-1)]
output_path = args[[length(args)]]

renv::load(project = code_directory)

library(tidyverse)
library(Seurat)
library(tidyseurat)
library(glue)

# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE)

input_paths |>
  map_dfr(readRDS) |>
  pull(variable_features) |>
  unique() |>
  saveRDS(output_path)
