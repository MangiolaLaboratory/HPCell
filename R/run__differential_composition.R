
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
library(patchwork)

load(glue("{code_directory}/data/theme_multipanel.rda"))

# Create dir
output_path_results |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

cell_data =
  input_annotation_label_transfer |>
  map_dfr(
    ~ readRDS(.x) |>
      mutate(file = basename(.x)) |>
      extract(file, "sample",  "(.+)__.+")
  ) |>
  left_join(readRDS(input_metadata), by="sample") |>

  # Transform time
  mutate(days_since_symptom_onset_log = days_since_symptom_onset |> log1p()) |>
  mutate(severity = as.character(severity)) |>  mutate(severity = if_else(severity =="NA", "healthy", severity)) |>

  # Filter out erhythrocites
  filter(!predicted.celltype.l2 %in% c("Eryth", "Platelet"))

estimates =
  cell_data |>

  # This will not be needed with contrasts
  sccomp_glm(
    formula_composition = ~ 0 + severity,
    formula_variability = ~ 0 + severity,
    .sample = sample,
    .cell_group = predicted.celltype.l2,
    contrasts = c("(1/2 * severitymoderate + 1/2 * severitysevere) - severityhealthy", "severitysevere - (1/2 * severitymoderate + 1/2 * severityhealthy)"),
    check_outliers = TRUE,
    bimodal_mean_variability_association = TRUE
  )


plots = plot_summary(estimates)


estimates_time =
  cell_data |>
  filter(!is.na(days_since_symptom_onset_log)) |>
  sccomp_glm(
    formula_composition = ~ days_since_symptom_onset_log + severity + donor,
    formula_variability = ~ days_since_symptom_onset_log + severity,
    .sample = sample,
    .cell_group = predicted.celltype.l2,
    check_outliers = TRUE,
    bimodal_mean_variability_association = TRUE
  )

plots_CI_time =
  estimates_time |>
  filter(parameter =="days_since_symptom_onset_log") |>

  # Reshape
  pivot_longer(c(contains("c_"), contains("v_")),names_sep = "_" , names_to=c("which", "estimate") ) |>
  #drop_na() |>
  pivot_wider(names_from = estimate, values_from = value) |>

  nest(data = -c(parameter, which)) |>
  mutate(plot = pmap(
    list(data, which, parameter),
    ~  ggplot(..1, aes(x=effect, y=fct_reorder(predicted.celltype.l2, effect))) +
      geom_vline(xintercept = 0.2, colour="grey") +
      geom_vline(xintercept = -0.2, colour="grey") +
      geom_errorbar(aes(xmin=lower, xmax=upper, color=FDR<0.025)) +
      geom_point() +
      scale_color_brewer(palette = "Set1") +
      xlab("Credible interval of the slope") +
      ylab("Cell group") +
      ggtitle(sprintf("%s %s", ..2, ..3)) +
      theme(legend.position = "bottom") +
      theme_bw()
  )) %>%
  pull(plot) |>
  wrap_plots(ncol=2) +

  # Style
  plot_layout(guides = 'collect') &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"), legend.position = "bottom")

# Plot
list(
  plots,
  plots_CI_time
) |>
  saveRDS(output_path_plot_credible_intervals)

# This will not be needed anymore
(
  plots$boxplot[[1]] &
    theme_multipanel
) |>
  saveRDS(output_path_plot_boxplot)

list(estimates, estimates_time) |> saveRDS(output_path_results)
