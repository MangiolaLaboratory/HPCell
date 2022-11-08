
set.seed(42)

number_features_overall = 300
number_features_per_cell_type = 300

# Read arguments
args = commandArgs()
code_directory = args[[1]]
input_paths = args[2:(length(args)-2)]
cell_type_column_for_subsetting = args[[length(args)-1]]
output_path = args[[length(args)]]

#renv::load(project = code_directory)

library(dplyr); library(tidyr); library(ggplot2)
library(Seurat)
library(tidyseurat)
library(glue)

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
      filter(!predicted.celltype.l2 %in% c("Eryth", "Platelet"))

  )) |>

  pull(seurat) |>
  reduce(merge) |>

  # Sample up to a plateau to avoid extreme cell_type bias
  nest(data = -sample) %>%
  mutate(n = map_int(data, ~ ncol(.x))) %>%
  mutate(upper_quantile = quantile(n, 0.75) %>% as.integer()) %>%
  mutate(data = map2(
    data, upper_quantile,
    ~ sample_n(.x, min(ncol(.x), .y), replace = FALSE)
  )) %>%
  filter(n>1) %>%
  unnest(data) |>

  # Normalise before - https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#normalizing-the-data-1
  NormalizeData(assay="RNA")


bind_rows(

  # General variable genes
  counts |>
  FindVariableFeatures(nfeatures = number_features_overall, assay="RNA") |>
  VariableFeatures(assay="RNA") |>
  as_tibble() |>
  rename("variable_features" = "value") |>
  mutate(group= "overall"),

# Per cell type
  counts |>

  # If cell_type_column_for_subsetting == null calculate clusters
  when(
    cell_type_column_for_subsetting=="none" ~
      counts |>
      FindVariableFeatures(nfeatures = number_features_overall, assay="RNA")  |>
      NormalizeData(assay="RNA") |>
      ScaleData(assay="RNA") |>
      RunPCA(assay="RNA") |>
      FindNeighbors(dims = 1:20) |>
      FindClusters(resolution = 0.5) |>
      nest(data = -seurat_clusters),

    # If I have label
    ~ nest(., data = -!!as.symbol(cell_type_column_for_subsetting))
  ) |>

  # Get feature within each cluster/cell-type
  mutate(variable_features = map(
    data,
    ~ .x |>
      FindVariableFeatures(nfeatures = number_features_per_cell_type, assay="RNA") |>
      VariableFeatures(assay="RNA")
  )) |>
  select(-data) |>
  unnest(variable_features) |>

  # Rename group column with cluster
  when(
    cell_type_column_for_subsetting!="none" ~ rename(., group := !!as.symbol(cell_type_column_for_subsetting) ),
    ~ rename(., group := seurat_clusters )
  )

) %>%
  saveRDS(output_path)
