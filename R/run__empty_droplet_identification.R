# Rscript analysis/split_single_cell_object_into_files.R data/alive_identification/counts_alive.rds data/alive_identification _alive

set.seed(42)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_path = args[[2]]
output_path_result = args[[3]]
output_path_plot_pdf = args[[4]]
output_path_plot_rds = args[[5]]

renv::activate(project = code_directory)

library(tidyverse)
library(purrr)
library(Seurat)
library(tidyseurat)
library(glue)
library(scater)
library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(here)



# Create dir
output_path_result |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

input_file = readRDS(input_path)


# Calculate bar-codes ranks
barcode_ranks = barcodeRanks(input_file)

barcode_table =
  barcode_ranks |>
  as_tibble(rownames = ".cell") |>
  mutate(
    knee =  metadata(barcode_ranks)$knee,
    inflection =  metadata(barcode_ranks)$inflection
  ) |>
  arrange(rank) |>
  mutate(empty_droplet = total<inflection)

barcode_table |>  saveRDS(output_path_result)

# Plot bar-codes ranks
plot_barcode_ranks =
  barcode_table %>%
  ggplot2::ggplot(aes(rank, total)) +
  geom_point() +
  geom_line(aes(rank, fitted), color="red") +
  geom_hline(aes(yintercept = knee), color="dodgerblue") +
  geom_hline(aes(yintercept = inflection), color="forestgreen") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()

plot_barcode_ranks |> saveRDS(output_path_plot_rds)

ggsave(
  output_path_plot_pdf,
  plot = plot_barcode_ranks,
  useDingbats=FALSE,
  units = c("mm"),
  width = 183/2 ,
  height = 183/2,
  limitsize = FALSE
)


