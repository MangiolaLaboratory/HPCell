
# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_pseudobulk = args[[2]]
input_metadata = args[[3]]
output_path_results = args[[4]]
output_path_plot_densities = args[[5]]
output_path_plot_significant = args[[6]]

# renv::load(project = code_directory)

library(dplyr); library(tidyr); library(ggplot2)
library(glue)
library(tidybulk)
library(tidySummarizedExperiment)
library(tidyHeatmap)

load(glue("{code_directory}/data/theme_multipanel.rda"))

# Create dir
output_path_results |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)



scaled_counts =
  readRDS(input_pseudobulk) |>

  left_join(
    readRDS(input_metadata)  |>
      mutate(severity = if_else(severity=="NA", "healthy", as.character(severity))) ,
    by="sample"
  ) |>

  identify_abundant(factor_of_interest = c(severity, predicted.celltype.l2)) |>
  scale_abundance(method="TMMwsp")


# Plot densities
plot_densities =
  scaled_counts |>
  filter(.abundant) |>

  # Calculate library size
  left_join(
    scaled_counts |> group_by(.sample) |>  summarise(count_sum = sum(count)),
    by = ".sample"
  ) |>

  ggplot(aes(count_scaled +1, group=sample, color=count_sum)) +
  geom_density() +
  scale_x_log10() +
  scale_color_gradient(trans = "log") +
  facet_wrap(~predicted.celltype.l2) +
  theme_multipanel

plot_densities |> saveRDS(output_path_plot_densities)


counts_tested =
  scaled_counts |>
  nest(se = -predicted.celltype.l2) |>

  # Filter if not enough samples per condition
  filter(map_int(
    se,
    ~ .x |>
      distinct(sample, severity) |>
      count(severity) |>
      pull(n) |>
      max()
  )>1) |>

  # test
  mutate(se = map(
    se,
    ~ .x %>% {print(.x); .x} |>
      test_differential_abundance(
        ~ 0 + severity ,
        contrasts = c("severitymoderate-severityNA", "severitysevere-severitymoderate"),
        method = "edger_robust_likelihood_ratio",
        test_above_log2_fold_change = 1,
        scaling_method = "TMMwsp"
      )
  ))

counts_tested |> saveRDS(output_path_results)


plot_significant_boxplot =

  counts_tested |>

  # Filter signifant
  mutate(significant = map(
    se,
    ~ .x |>
      filter(.abundant) |>
      filter(`FDR___severitymoderate-severityNA` < 0.05 | `FDR___severitysevere-severitymoderate` < 0.05) |>
      pivot_longer(names_sep = "___", cols = contains("___"), names_to = c("statistics", "contrast"), values_to = "value")
  )) |>
  select(-se) |>
  unnest(significant) |>

  # Plot
  ggplot(aes(.feature,count_scaled + 1, fill=fct_relevel(severity, c("NA", "moderate", "severe")))) +
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(~predicted.celltype.l2, scale="free_x", ncol=1) +
  theme_multipanel +
  theme(axis.text.x = element_text(angle=90))


plot_significant_heatmap =
  counts_tested |>

  # Filter signifant
  mutate(significant = map(
    se,
    ~ .x |>
      filter(.abundant) |>
      filter(`FDR___severitymoderate-severityNA` < 0.05 | `FDR___severitysevere-severitymoderate` < 0.05) |>
      as_tibble()
  )) |>
  select(-se) |>
  unnest(significant)  |>
  tidybulk:::drop_class("tidySummarizedExperiment_nested") |>
  unite(".feature", c(.feature, predicted.celltype.l2), remove = FALSE)  |>
  extract(.sample, ".sample", "(.+)___.+") |> group_by(predicted.celltype.l2) |>
  heatmap(.feature, .sample, count_scaled, transform = log1p) |>
  add_tile(severity)

list(plot_significant_boxplot, plot_significant_heatmap) |> saveRDS(output_path_plot_significant)
