set.seed(42)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_paths = args[2:(length(args)-1)]
output_path = args[[length(args)]]

renv::load(project = code_directory)

library(dplyr); library(tidyr); library(ggplot2)
library(Seurat)
library(tidyseurat)
library(glue)
library(purrr)
library(tibble)

# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE)

number_of_files = length(input_paths)

variable_features_df =
  input_paths |>
  enframe(value = "file") |>
  select(-name) |>
  mutate(data = map(file, readRDS)) |>
  unnest(data)

gene_intersection =
  variable_features_df |>
  filter(group == "all_features") |>
  count(feature) |>
  filter(n == number_of_files) |>
  pull(feature)

# ERROR
stopifnot(
  "There is 0 features that is common to all samples wirthin the RNA assay. This should not happen." =
    length(gene_intersection) > 0

  )

# Select only the genes present in all
variable_features_df =
  variable_features_df |>
  filter(feature %in% gene_intersection)


subset_top_rank_variable_genes_across_batches(
   variable_features_df |> filter(group == "variable_overall"),
   variable_features_df |> filter(!group %in% c("variable_overall", "all_features")),
   group,
   file,
    features_number_independent_of_cell_groups = 2000,
    features_number_per_cell_group = 300
) |>
  saveRDS(output_path)


