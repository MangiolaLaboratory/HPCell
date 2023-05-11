# Rscript analysis/split_single_cell_object_into_files.R data/alive_identification/counts_alive.rds data/alive_identification _alive

set.seed(42)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_files = args[2:(length(args)-3)]
input_path_10x_raw = args[2:(length(args)-3)]
output_path_result = args[[length(args)-2]]
output_path_plot_pdf = args[[length(args)-1]]
output_path_plot_rds = args[[length(args)]]

input_path_10x_raw = input_files[1:(length(input_files)/2)]
input_path_demultiplexed = args[[(length(input_files)/2+1):length(input_files)]]

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
library(celda)
library(scater)
library(tidySingleCellExperiment)
library(DropletUtils)


# Create dir
output_path_result |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

"~/PostDoc/covid19pbmc/data/C120_COVID_PBMC_batch14/extdata/CellRanger/C120_batch14_1/outs/multi/count/raw_feature_bc_matrix/"

raw = input_path_10x_raw |> map(read10xCounts) |> reduce(bind_rows)
filtered = input_path_demultiplexed |> map(~ readRDS(.x) |> mutate(file = .x)) |> reduce(bind_rows)

colnames(raw) = colData(raw)$Barcode
colnames(filtered) = colData(filtered)$Barcode
raw = raw[,!colnames(raw) %in% colnames(filtered)]

# Drop top 10% of background
raw[,colSums(raw)<quantile(colSums(raw), 0.95)]


filtered <-
  decontX(x = filtered, background = raw) |>
  as.Seurat(counts = "decontXcounts", data = NULL) |>
  nest(data = -file) |>
  mutate(saved = map2(data, file, ~ .x |> saveRDS(glue("DIR_OUT{.y}"))))

# umap <- reducedDim(pbmc4k, "decontX_UMAP")
# plotDimReduceCluster(x = pbmc4k$,
#                      dim1 = umap[, 1], dim2 = umap[, 2])
# plotDecontXContamination(pbmc4k)
#
# pbmc4k <- logNormCounts(pbmc4k)
# plotDimReduceFeature(as.matrix(logcounts(pbmc4k)),
#                      dim1 = umap[, 1],
#                      dim2 = umap[, 2],
#                      features = c("CD3D", "CD3E", "GNLY",
#                                   "LYZ", "S100A8", "S100A9",
#                                   "CD79A", "CD79B", "MS4A1"),
#                      exactMatch = TRUE)



