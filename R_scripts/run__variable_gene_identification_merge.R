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

variable_within_cell_types =
  variable_features_df |>
  filter(!group %in% c("variable_overall", "all_features")) |>
  count(feature, group) |>
  with_groups(group, ~ .x |> arrange(desc(n)) |> slice(1:300))

variable_across_cell_types =
  variable_features_df |>

  # Filter files that have more than 10 cell types + variable_overall
  nest(data = -file) |>
  filter(map_int(data, ~ .x |> distinct(group) |> nrow()) > 11) |>
  unnest(data) |>

  filter(group == "variable_overall") |>
  count(feature, group) |>
  with_groups(group, ~ .x |> arrange(desc(n)) |> slice(1:2000))


bind_rows(

    # Cell type specific
    variable_within_cell_types,

    # Overall
    variable_across_cell_types

  ) |>
  pull(feature) |>
  unique() |>
  saveRDS(output_path)
