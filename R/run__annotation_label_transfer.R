# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

# Rscript /home/users/allstaff/mangiola.s/PostDoc/jascap/R/run__annotation_label_transfer.R
# input_path_demultiplexed = "data/3_prime_batch_1/input_files/0483-002_NA.rds"
# input_path_empty_droplets = "data/3_prime_batch_1/preprocessing_results/empty_droplet_identification/0483-002_NA__empty_droplet_identification_output.rds"
# reference_azimuth_path = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/reference_azimuth.rds"
# output_path = "data/3_prime_batch_1/preprocessing_results/annotation_label_transfer/0483-002_NA__annotation_label_transfer_output.rds"


# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_path_demultiplexed = args[[2]]
input_path_empty_droplets = args[[3]]
reference_azimuth_path = args[[4]]
output_path = args[[5]]

renv::load(project = code_directory)

library(dplyr); library(tidyr); library(ggplot2)
library(Seurat)
library(tidyseurat)
library(glue)
library(SingleR)
library(celldex)
library(scuttle)
library(purrr)

# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

# Load reference PBMC
# reference_azimuth <- LoadH5Seurat("data//pbmc_multimodal.h5seurat")
# reference_azimuth |> saveRDS("analysis/annotation_label_transfer/reference_azimuth.rds")

reference_azimuth = readRDS(reference_azimuth_path)

# Reading input
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
azimuth_annotation =
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
   select(.cell, predicted.celltype.l1, predicted.celltype.l2, contains("refUMAP"))

rm(input_file)
gc()

# SingleR




  input_file_sce =
  readRDS(input_path_demultiplexed) |>

  # Filter empty
  left_join(readRDS(input_path_empty_droplets), by = ".cell") |>
  tidyseurat::filter(!empty_droplet) |>
  as.SingleCellExperiment() |>
    logNormCounts()

  # Blueprint
  blueprint <- BlueprintEncodeData()
  blueprint_annotation_fine =
    input_file_sce |>
   SingleR(ref = blueprint,
           assay.type.test=1,
           labels = blueprint$label.fine
          ) |>
    as_tibble(rownames = ".cell") |>
    select(.cell, blueprint_first.labels.fine = first.labels)

  blueprint_annotation_coarse =
    input_file_sce |>
    SingleR(ref = blueprint,
            assay.type.test=1,
            labels = blueprint$label.main
    ) |>
    as_tibble(rownames = ".cell") |>
    select(.cell, blueprint_first.labels.coarse = first.labels)

  rm(blueprint)
  gc()

  # Monaco
  MonacoImmuneData = MonacoImmuneData()
  monaco_annotation_fine =
    input_file_sce |>
    SingleR(ref = MonacoImmuneData,
            assay.type.test=1,
            labels = MonacoImmuneData$label.fine
          ) |>
    as_tibble(rownames = ".cell") |>
    select(.cell, monaco_first.labels.fine = first.labels)

  monaco_annotation_coarse =
    input_file_sce |>
    SingleR(ref = MonacoImmuneData,
            assay.type.test=1,
            labels = MonacoImmuneData$label.main
    ) |>
    as_tibble(rownames = ".cell") |>
    select(.cell, monaco_first.labels.coarse = first.labels)

  rm(MonacoImmuneData)
  gc()

  # Clean
  rm(input_file_sce)
  gc()

  # Join annotation
  azimuth_annotation |>
    left_join(blueprint_annotation_fine) |>
    left_join(blueprint_annotation_coarse) |>
    left_join(monaco_annotation_fine) |>
    left_join(monaco_annotation_coarse) |>

   # Save
   saveRDS(output_path)
