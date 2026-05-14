library(testthat)
library(HPCell)
library(Seurat)
library(scRNAseq)
library(magrittr)

tissue <- "pbmc"
RNA_assay_name <- "originalexp"

input_seurat_abc <-
  HeOrganAtlasData(ensembl = FALSE, location = FALSE) |>
  as.Seurat(data = NULL) |>
  subset(subset = Tissue %in% c("Blood"))

reference_label_fine <- HPCell:::reference_label_fine_id(tissue)

empty_droplets_tbl <-
  HPCell:::empty_droplet_id(
    input_seurat_abc,
    feature_nomenclature = "symbol",
    filter_empty_droplets = TRUE
  )

annotation_label_transfer_tbl <-
  HPCell:::annotation_label_transfer(input_seurat_abc, empty_droplets_tbl)

alive_identification_tbl <-
  HPCell:::alive_identification(
    input_seurat_abc,
    empty_droplets_tbl,
    annotation_label_transfer_tbl
  )

doublet_identification_tbl <-
  HPCell:::doublet_identification(
    input_seurat_abc,
    empty_droplets_tbl,
    alive_identification_tbl
  )

cell_cycle_score_tbl <-
  HPCell:::cell_cycle_scoring(input_seurat_abc, empty_droplets_tbl)

non_batch_variation_removal_S <-
  HPCell:::non_batch_variation_removal(
    input_seurat_abc,
    empty_droplets_tbl,
    alive_identification_tbl,
    cell_cycle_score_tbl,
    assay = NULL
  )

preprocessing_output_S <-
  HPCell:::preprocessing_output(
    input_seurat_abc,
    empty_droplets_tbl,
    non_batch_variation_removal_S,
    alive_identification_tbl,
    cell_cycle_score_tbl,
    annotation_label_transfer_tbl,
    doublet_identification_tbl
  )

test_that("input_seurat_works", {
  expect_s4_class(input_seurat_abc, "Seurat")
})

test_that("reference_label_fine works", {
  expect_true(
    "monaco_first.labels.fine" %in% reference_label_fine ||
      "blueprint_first.labels.fine" %in% reference_label_fine ||
      "none" %in% reference_label_fine
  )
})

test_that("empty_droplets_works", {
  expect_s3_class(empty_droplets_tbl, "tbl_df")
  expect_lte(nrow(empty_droplets_tbl), ncol(input_seurat_abc))
})

test_that("cell_cycle_score_works", {
  expect_s3_class(cell_cycle_score_tbl, "tbl_df")
  expected_colnames <- c("S.Score", "G2M.Score", "Phase")
  expect_true(all(expected_colnames %in% colnames(cell_cycle_score_tbl)))
})

test_that("annotation_label_transfer_works", {
  expect_equal(ncol(annotation_label_transfer_tbl), 9)
})

test_that("alive_identification_works", {
  expect_s3_class(alive_identification_tbl, "tbl_df")
  expected_colnames <- c(
    "subsets_Mito_sum",
    "subsets_Mito_detected",
    "subsets_Mito_percent"
  )
  expect_true(all(expected_colnames %in% colnames(alive_identification_tbl)))
})

test_that("non_batch_variation_removal_S_dimensions", {
  num_features_input <- nrow(input_seurat_abc)
  num_cells_input <- ncol(input_seurat_abc)

  num_features_non_batch <- nrow(non_batch_variation_removal_S)
  num_cells_non_batch <- ncol(non_batch_variation_removal_S)

  expect_true(num_features_non_batch <= num_features_input)
  expect_true(num_cells_non_batch <= num_cells_input)
})

test_that("Doublet_identification_works", {
  expect_s3_class(doublet_identification_tbl, "tbl_df")
  expected_colnames <- c(
    "scDblFinder.class",
    "scDblFinder.score",
    "scDblFinder.weighted",
    "scDblFinder.cxds_score"
  )
  expect_true(all(expected_colnames %in% colnames(doublet_identification_tbl)))
})

test_that("Preprocessing_works", {
  expect_s4_class(preprocessing_output_S, "Seurat")
})
