
# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_path_preprocessing = args[[2]]
output_path = args[[3]]

renv::activate(project = code_directory)

library(tidyverse)
library(Seurat)
library(glue)
library(CellChat)
library(tidyseurat)
library(tidySingleCellExperiment)
library(furrr)
library(stringr)

plan(multisession, workers=3)

# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

counts_cellchat =
  input_path_preprocessing |>
  readRDS() |>
  createCellChat(group.by = "predicted.celltype.l2") |>
  setIdent( ident.use = "predicted.celltype.l2")


tibble(DB = c("Secreted Signaling", "ECM-Receptor" , "Cell-Cell Contact" )) |>
  mutate(sample =  input_path_preprocessing |> basename() |> str_remove("__.+")) |>
  mutate(data = list(counts_cellchat)) |>
  mutate(data = future_map2(
    data, DB,
    ~ {
      print(.y)
      .x@DB <- subsetDB(CellChatDB.human, search = .y)

      .x |>
        subsetData() |>
        identifyOverExpressedGenes() |>
        identifyOverExpressedInteractions() |>
        projectData(PPI.human) |>
        computeCommunProb() |>
        filterCommunication() |>
        computeCommunProbPathway() |>
        aggregateNet()

    }
  )) |>

  # Save
  saveRDS(output_path)

