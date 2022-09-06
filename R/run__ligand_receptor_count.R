
# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_path_preprocessing = args[[2]]
output_path = args[[3]]

renv::load(project = code_directory)

library(tidyverse)
library(Seurat)
library(glue)
library(CellChat)
library(tidyseurat)
library(tidySingleCellExperiment)
library(stringr)

#library(furrr)
#plan(multisession, workers=3)

# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

counts_cellchat =
  input_path_preprocessing |>
  readRDS() |>
  createCellChat(group.by = "predicted.celltype.l2") |>
  setIdent( ident.use = "predicted.celltype.l2")

# More robust implementation that does not fail if no results
computeCommunProbPathway = function (object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05)
{
  if (is.null(net)) {
    net <- object@net
  }
  if (is.null(pairLR.use)) {
    pairLR.use <- object@LR$LRsig
  }
  prob <- net$prob
  prob[net$pval > thresh] <- 0
  pathways <- unique(pairLR.use$pathway_name)
  group <- factor(pairLR.use$pathway_name, levels = pathways)
  prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum),
                         c(2, 3, 1))
  pathways.sig <- pathways[apply(prob.pathways, 3, sum) !=
                             0]
  prob.pathways.sig <- prob.pathways[, , pathways.sig, drop=FALSE]
  idx <- sort(apply(prob.pathways.sig, 3, sum), decreasing = TRUE,
              index.return = TRUE)$ix
  pathways.sig <- pathways.sig[idx]
  prob.pathways.sig <- prob.pathways.sig[, , idx]
  if (is.null(object)) {
    netP = list(pathways = pathways.sig, prob = prob.pathways.sig)
    return(netP)
  }
  else {
    object@netP$pathways <- pathways.sig
    object@netP$prob <- prob.pathways.sig
    return(object)
  }
}

tibble(DB = c("Secreted Signaling", "ECM-Receptor" , "Cell-Cell Contact" )) |>
  mutate(sample =  input_path_preprocessing |> basename() |> str_remove("__.+")) |>
  mutate(data = list(counts_cellchat)) |>
  mutate(data = map2(
    data, DB,
    ~ {
      print(.y)
      .x@DB <- subsetDB(CellChatDB.human, search = .y)

      x = .x |>
        subsetData() |>
        identifyOverExpressedGenes() |>
        identifyOverExpressedInteractions() |>
        projectData(PPI.human)
browser()
    if(nrow(x@LR$LRsig)==0) return(NA)

     x |>
        computeCommunProb(do.fast = FALSE) |>
        filterCommunication() |>
        computeCommunProbPathway() |>
        aggregateNet()

    }
  )) |>

  # Save
  saveRDS(output_path)

