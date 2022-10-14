# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_metadata = args[[2]]
input_files = args[3:(length(args)-2)]
output_path_dataframe = args[[length(args)-1]]
output_path_plot_umap = args[[length(args)]]

# Divide input
input_demulatiplexed = input_files[1:(length(input_files)/2)]
input_empty_droplets = input_files[((length(input_files)/2)+1):length(input_files)]

renv::load(project = code_directory)

library(tidyverse)
library(Seurat)
library(tidyseurat)
library(glue)
library(purrr)
library(ggupset)
load(glue("{code_directory}/data/theme_multipanel.rda"))

assay_of_choice = "RNA"

# Create dir
output_path_dataframe |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

metadata =

variable_genes_per_sample =

  # Reading input
  tibble(
    seurat_obj_list = input_demulatiplexed |> map(readRDS),
    empty_droplets_obj_list = input_empty_droplets |> map(readRDS)
  ) |>

  # Filter
  mutate(seurat_obj_list = map2(seurat_obj_list, empty_droplets_obj_list, ~ .x |> left_join(.y) |> filter(!empty_droplet))) |>
  select(seurat_obj_list) |>
  tidyseurat:::unnest.tidyseurat_nested(seurat_obj_list) |>
  ScaleData(assay=assay_of_choice, return.only.var.genes=FALSE)  |>

  # Add metadata
  left_join(readRDS(input_metadata)) |>

  # Get variable genes
  nest(data = -sample) |>
  mutate(variable_genes = map(
    data,
    ~ .x |>
      FindVariableFeatures(assay=assay_of_choice, nfeatures = 500) |>
      VariableFeatures(assay=assay_of_choice)
  ))  |>

  # Subset
  mutate(n = map_int(data, ~ ncol(.x))) %>%
  mutate(upper_quantile = quantile(n, 0.50) %>% as.integer()) %>%
  mutate(data = map2(
    data, upper_quantile,
    ~ sample_n(.x, min(ncol(.x), .y), replace = FALSE)
  ))



data_umap =
  variable_genes_per_sample |>
  select(sample, data) |>
  unnest(data) %>%
  ScaleData(assay=assay_of_choice, return.only.var.genes=FALSE) %>%
  # Variable genes
  {
    .x = (.)
    VariableFeatures(.x) =
      variable_genes_per_sample |>
      select(sample, variable_genes) |>
      unnest(variable_genes) |>
      add_count(variable_genes) |>
      filter(n>3) |>
      pull(variable_genes) |>
      unique()
    .x
  } |>

  # UMAP
  RunPCA(dims = 1:30, assay=assay_of_choice) |>
  RunUMAP(dims = 1:30, spread = 0.5,min.dist  = 0.01, n.neighbors = 10L) |>
  as_tibble()

# Plot
plot_sample_color =
  data_umap |>

  # I HAVE TO DROP THIS
  mutate(batch = 1) |>

  # UMAP
  ggplot(aes(UMAP_1, UMAP_2, color = sample)) +
  geom_point(size=0.2) +
  facet_wrap(~batch) +
  guides(color="none") +
  theme_multipanel

plot_severity_color =
  data_umap |>

  # I HAVE TO DROP THIS
  mutate(batch = 1) |>

  # UMAP
  ggplot(aes(UMAP_1, UMAP_2, color = severity)) +
  geom_point(size=0.2) +
  facet_wrap(~batch) +
  guides(color="none") +
  theme_multipanel

plot_batch_color =
  data_umap |>

  # I HAVE TO DROP THIS
  mutate(batch = 1) |>

  # UMAP
  ggplot(aes(UMAP_1, UMAP_2, color = batch)) +
  geom_point(size=0.2) +
  theme_multipanel


 # Save
saveRDS(data_umap, output_path_dataframe)
saveRDS(list(
  plot_sample_color,
   plot_severity_color,
  plot_batch_color
), output_path_plot_umap) +
  theme_multipanel
