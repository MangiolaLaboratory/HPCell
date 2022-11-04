# Rscript analysis/split_single_cell_object_into_files.R data/alive_identification/counts_alive.rds data/alive_identification _alive

set.seed(42)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_path = args[[2]]
output_path_result = args[[3]]
output_path_plot_pdf = args[[4]]
output_path_plot_rds = args[[5]]

renv::load(project = code_directory)

library(dplyr); library(tidyr); library(ggplot2)
library(purrr)
library(Seurat)
library(tidyseurat)
library(glue)
library(scater)
library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(here)
library(stringr)

significance_threshold = 0.001

# Create dir
output_path_result |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

input_file = readRDS(input_path)

# Genes to exclude
location <- mapIds(
  EnsDb.Hsapiens.v86,
  keys=rownames(input_file),
  column="SEQNAME",
  keytype="SYMBOL"
)
mitochondrial_genes = which(location=="MT") |> names()
ribosome_genes = rownames(input_file) |> str_subset("^RPS|^RPL")

# Calculate bar-codes ranks
barcode_ranks = barcodeRanks(input_file)

# Set the minimum total RNA per cell for ambient RNA
if(min(barcode_ranks$total) < 100) { lower = 100 } else {
  lower = quantile(barcode_ranks$total, 0.05)

  write_lines(
    glue("{input_path} has supposely empty droplets with a lot of RNAm maybe a lot of ambient RNA? Please investigate"),
    file = glue("{dirname(output_path_result)}/warnings_emptyDrops.txt"),
    append = T
  )
}


# Remove genes from input
barcode_table =

  # Proper classification
  input_file@assays$RNA@counts[!rownames(input_file@assays$RNA@counts) %in% c(mitochondrial_genes, ribosome_genes),, drop=FALSE] |>
  emptyDrops( test.ambient = TRUE, lower=lower) |>
  as_tibble(rownames = ".cell") |>
  mutate(empty_droplet = !FDR< significance_threshold) |>
  replace_na(list(empty_droplet = TRUE)) |>

  # barcode ranks
  left_join(
    barcode_ranks |>
      as_tibble(rownames = ".cell") |>
      mutate(
        knee =  metadata(barcode_ranks)$knee,
        inflection =  metadata(barcode_ranks)$inflection
      )
  )


barcode_table |>  saveRDS(output_path_result)

# Plot bar-codes ranks
plot_barcode_ranks =
  barcode_table %>%
  ggplot2::ggplot(aes(rank, total)) +
  geom_point(aes(color = empty_droplet, size = empty_droplet )) +
  geom_line(aes(rank, fitted), color="purple") +
  geom_hline(aes(yintercept = knee), color="dodgerblue") +
  geom_hline(aes(yintercept = inflection), color="forestgreen") +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
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


