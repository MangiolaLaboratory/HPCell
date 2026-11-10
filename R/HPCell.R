#' HPCell: Massively-Parallel R Native Pipeline for Single-Cell Analysis
#'
#' @importFrom methods is
#' @importFrom stats formula quantile setNames
#' @importFrom utils data head tail
"_PACKAGE"

# Register tidy* S3 methods (e.g. left_join.Seurat, dplyr verbs on SCE)
# without import(tidy*) in NAMESPACE (avoids export conflicts with dplyr/broom on check).
.onLoad <- function(libname, pkgname) {
  loadNamespace("tidyseurat")
  loadNamespace("tidySingleCellExperiment")
  loadNamespace("tidySummarizedExperiment")
}

.myDataEnv <- new.env(parent = emptyenv()) # not exported

.data_internal <- function(dataset) {
  if (!exists(dataset, envir = .myDataEnv)) {
    utils::data(list = c(dataset), envir = .myDataEnv)
  }
}