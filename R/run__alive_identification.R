# Rscript analysis/split_single_cell_object_into_files.R data/alive_identification/counts_alive.rds data/alive_identification _alive

# This script slipts the dataset and creates the list of files in a specific directory

set.seed(42)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_path_demultiplexed = args[[2]]
input_path_empty_droplets = args[[3]]
output_path = args[[4]]

renv::load(project = code_directory)

library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(purrr)
library(Seurat)
library(tidyseurat)
library(glue)
library(scater)



# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

input_file =
  readRDS(input_path_demultiplexed) |>
  left_join(readRDS(input_path_empty_droplets), by=".cell") |>
  tidyseurat::filter(!empty_droplet)

location <- mapIds(
  EnsDb.Hsapiens.v86,
  keys=rownames(input_file),
  column="SEQNAME",
  keytype="SYMBOL"
)

mitochondrion =
  input_file |>
  GetAssayData( slot = "counts", assay="RNA") |>

  # Join mitochondrion statistics
  scuttle::perCellQCMetrics(subsets=list(Mito=which(location=="MT"))) |>
  as_tibble(rownames = ".cell") |>
  dplyr::select(-sum, detected) |>

  # Label cells
  mutate(high_mitochondrion = isOutlier(subsets_Mito_percent, type="higher"))

ribosome =

  input_file |>

  tidyseurat::select(.cell) |>

  # Join mitochondrion statistics
  mutate(mito_RPS = PercentageFeatureSet(input_file,  pattern = "^RPS|^RPL", assay = "RNA")[,1]  ) |>

  # Label cells
  mutate(high_RPS = isOutlier(mito_RPS, type="higher")) |>

  as_tibble() |>
  dplyr::select(.cell, mito_RPS, high_RPS)

# Save
mitochondrion |>
  left_join(ribosome, by=".cell")|>
  saveRDS(output_path)

# plot_QC =
#   input_file |>
#   left_join(annotation, by=".cell") |>
#   pivot_longer(c(nFeature_RNA , nCount_RNA, subsets_Mito_percent, mito_RPS)) |>
#   ggplot(aes(x=sample,y=value,
#              color = high_mitochondrion,
#              alpha=high_mitochondrion,
#              size= high_mitochondrion
#   )) +
#   ggbeeswarm::geom_quasirandom() +
#   facet_wrap(~name, scale="free_y") +
#
#   # Customisation
#   scale_color_manual(values=c("black", "#e11f28")) +
#   scale_size_discrete(range = c(0, 2)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
#
# ggsave(
#   "benign_mitochondrial.pdf",
#   plot = plot_QC,
#   useDingbats=FALSE,
#   units = c("mm"),
#   width = 183 ,
#   height = 120,
#   limitsize = FALSE
# )

