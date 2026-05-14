library(testthat)
library(HPCell)
library(Seurat)
library(scRNAseq)
library(crew)

test_that("HPCell pipeline works for the Seurat test data", {
  td <- file.path(
    tempdir(),
    paste0("hpcell_", format(Sys.time(), "%Y%m%d%H%M%S"))
  )
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  file <- file.path(td, "file.rds")
  scRNAseq::HeOrganAtlasData(
    tissue = "Blood",
    ensembl = FALSE,
    location = FALSE
  ) |>
    as.Seurat(data = NULL) |>
    saveRDS(file)

  store <- file.path(td, "_targets")

  c(my_sample = file) |>
    initialise_hpc(
      gene_nomenclature = "symbol",
      data_container_type = "seurat_rds",
      store = store,
      computing_resources = crew_controller_local(workers = 1)
    ) |>
    remove_empty_DropletUtils() |>
    annotate_cell_type() |>
    remove_dead_scuttle() |>
    score_cell_cycle_seurat() |>
    remove_doublets_scDblFinder() |>
    normalise_abundance_seurat_SCT(
      factors_to_regress = c(
        "subsets_Mito_percent",
        "subsets_Ribo_percent",
        "G2M.Score"
      )
    ) |>
    get_single_cell() |>
    evaluate_hpc()

  preprocessing_output <- targets::tar_read(single_cell, store = store)[[1]]

  metadata_columns <-
    if (inherits(preprocessing_output, "Seurat")) {
      colnames(preprocessing_output[[]])
    } else {
      colnames(SummarizedExperiment::colData(preprocessing_output))
    }

  metadata_columns |>
    stringr::str_which(pattern = "subsets_Mito_percent|subsets_Ribo_percent") |>
    expect_length(2)
})
