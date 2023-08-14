
set.seed(42)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
tissue = args[[2]]
input_paths = args[3:(length(args)-2)]
cell_type_column_for_subsetting = args[[length(args)-1]]
output_path = args[[length(args)]]

# renv::load(project = code_directory)

library(dplyr); library(tidyr); library(ggplot2)
library(Seurat)
library(tidyseurat)
library(glue)
library(purrr)
library(magrittr)
library(tibble)
library(jascap)
library(stringr)

# Divide into 5 chunks
input_paths_chunks = input_paths |> split( rep_len(1:5, length(input_paths)) |> sort())

# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE)

# Filter data
counts =
  tibble(
  input_data_paths = input_paths_chunks[[1]],
  input_empty_droplets_paths = input_paths_chunks[[2]],
  input_dead_cells_paths = input_paths_chunks[[3]],
  input_doublets_paths = input_paths_chunks[[4]],
  input_cell_type_annotation_paths = input_paths_chunks[[5]]
) |>
  mutate(seurat = pmap(
    list(input_data_paths, input_empty_droplets_paths, input_dead_cells_paths, input_doublets_paths, input_cell_type_annotation_paths),
    ~ readRDS(..1) |>

      # Filter empty droplets
      left_join(readRDS(..2) |> select(.cell, empty_droplet), by = ".cell") |>
      filter(!empty_droplet) |>

      # Filter dead cells
      left_join(readRDS(..3) |> select(.cell, alive), by = ".cell") |>
      filter(alive) |>

      # Filter doublets
      left_join(readRDS(..4) |> select(.cell, scDblFinder.class), by = ".cell") |>
      filter(scDblFinder.class=="singlet") |>

      # Filter Red blood cells and platelets using pbmc seurat for any tissue
      left_join(readRDS(..5), by = ".cell") |>

      # If PBMC filter for erythrocytes
      when(
        tissue |>
          tolower() |>
          equals("pbmc") ~
          filter(., !predicted.celltype.l2 %in% c("Eryth", "Platelet")),
        ~ (.)
      )


  )) |>

  pull(seurat) |>
  reduce(merge)


all_features_df =
  counts |>
  Assays("RNA") |>
  rownames() |>
  enframe(value = "feature") |>
  select(-name) |>
  mutate(group= "all_features")

# I BROKE cell_type_column_for_subsetting = "none"
seurat_to_variable_features(
    counts,
    "RNA",
    sample,
    sym(cell_type_column_for_subsetting),
    features_number_independent_of_cell_groups = 300,
    features_number_per_cell_group = 300
)|>
    bind_rows(all_features_df) |>
    saveRDS(output_path)


