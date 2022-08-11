# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_path_demultiplexed = args[[2]]
input_path_empty_droplets = args[[3]]
reference_azimuth_path = args[[4]]
output_path = args[[5]]

renv::load(project = code_directory)

library(tidyverse)
library(Seurat)
library(tidyseurat)
library(glue)

# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

# Load reference PBMC
# reference_azimuth <- LoadH5Seurat("data//pbmc_multimodal.h5seurat")
# reference_azimuth |> saveRDS("analysis/annotation_label_transfer/reference_azimuth.rds")

reference_azimuth = readRDS(reference_azimuth_path)

# Readingh input
input_file =

  readRDS(input_path_demultiplexed) |>

  # Filter empty
  left_join(readRDS(input_path_empty_droplets), by = ".cell") |>
  tidyseurat::filter(!empty_droplet) |>

  # Normalise RNA - not informed by smartly selected variable genes
  SCTransform(assay="RNA") |>

  # Normalise antibodies
  when(
    "ADT" %in% names(.@assays) ~ NormalizeData(., normalization.method = 'CLR', margin = 2, assay="ADT") ,
    ~ (.)
  )

# Define common anchors
anchors <- FindTransferAnchors(
  reference = reference_azimuth,
  query = input_file,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:30
)

# Mapping
 MapQuery(
  anchorset = anchors,
  query = input_file,
  reference = reference_azimuth ,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
) |>

   # Just select essential information
   as_tibble() |>
   select(.cell, predicted.celltype.l1, predicted.celltype.l2, contains("refUMAP")) |>

   # Save
   saveRDS(output_path)
