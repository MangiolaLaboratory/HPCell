
# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_metadata = args[[2]]
input_annotation_label_transfer = args[3:(length(args)-3)]

output_path_results = args[[length(args)-2]]
output_path_plot_credible_intervals = args[[length(args)-1]]
output_path_plot_boxplot = args[[length(args)]]


renv::load(project = code_directory)

library(tidyverse)
library(glue)
library(sccomp)
library(stringr)

load(glue("{code_directory}/data/theme_multipanel.rda"))

# Create dir
output_path_results |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

estimates1 =
  input_annotation_label_transfer |>
  map_dfr(
    ~ readRDS(.x) |>
      mutate(file = basename(.x)) |>
      extract(file, "sample",  "(.+)__.+")
    ) |>

  left_join(readRDS(input_metadata), by="sample") |>

  # This will not be needed with contrasts
  dplyr::filter(severity %in% c("NA", "moderate")) |>
  droplevels() |>
  mutate(severity = fct_relevel(severity, c("NA", "moderate"))) |>

  sccomp_glm(
    formula_composition = ~ severity ,
    formula_variability = ~ severity,
    .sample = sample,
    .cell_group = predicted.celltype.l2,
    check_outliers = TRUE,
    bimodal_mean_variability_association = TRUE
  )



plots1 = plot_summary(estimates1)



estimates2 =
  input_annotation_label_transfer |>
  map_dfr(
    ~ readRDS(.x) |>
      mutate(file = basename(.x)) |>
      extract(file, "sample",  "(.+)__.+")
  ) |>

  left_join(readRDS(input_metadata), by="sample") |>

  # This will not be needed with contrasts
  dplyr::filter(severity %in% c("severe", "moderate")) |>
  droplevels() |>
  mutate(severity = fct_relevel(severity, c("moderate", "severe"))) |>

  sccomp_glm(
    formula_composition = ~ severity ,
    formula_variability = ~ severity,
    .sample = sample,
    .cell_group = predicted.celltype.l2,
    check_outliers = TRUE,
    bimodal_mean_variability_association = TRUE
  )


plots2 = plot_summary(estimates2)

# Plot
(
  plots1$credible_intervals_1D  /
    plots2$credible_intervals_1D  &
    theme_multipanel
) |>
  saveRDS(output_path_plot_credible_intervals)

(
  plots1$boxplot[[1]]  +
    plots2$boxplot[[1]] &
    theme_multipanel
) |>
  saveRDS(output_path_plot_boxplot)

list(estimates1, estimates2) |> saveRDS(output_path_results)
