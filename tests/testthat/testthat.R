library(testthat)
library(HPCell)

## Define arguments 
filtered <- "TRUE"
tissue <- "pbmc"
RNA_assay_name<- "originalexp"
input_data_path<- "~/Documents/test_pipeline/file7de01ac8a860.rds"
input_file<- readRDS(input_data_path)
reference_azimuth = NULL
## Defining functions 
input_read_RNA_assay = add_RNA_assay(input_file, RNA_assay_name)
reference_label_fine = reference_label_fine_id(tissue)
empty_droplets_tbl = empty_droplet_id(input_read_RNA_assay, filtered)
annotation_label_transfer_tbl = annotation_label_transfer(input_read_RNA_assay,
                                                          reference_azimuth,
                                                          empty_droplets_tbl)
alive_identification_tbl = alive_identification(input_read_RNA_assay,
                                                empty_droplets_tbl,
                                                annotation_label_transfer_tbl)
doublet_identification_tbl = doublet_identification(input_read_RNA_assay,
                                                    empty_droplets_tbl,
                                                    alive_identification_tbl,
                                                    annotation_label_transfer_tbl,
                                                    reference_label_fine)
cell_cycle_score_tbl = cell_cycle_scoring(input_read_RNA_assay,
                                          empty_droplets_tbl)
non_batch_variation_removal_S = non_batch_variation_removal(input_read_RNA_assay,
                                                            empty_droplets_tbl,
                                                            alive_identification_tbl,
                                                            cell_cycle_score_tbl)
preprocessing_output_S = preprocessing_output(tissue,
                                              non_batch_variation_removal_S,
                                              alive_identification_tbl,
                                              cell_cycle_score_tbl,
                                              annotation_label_transfer_tbl,
                                              doublet_identification_tbl)
pseudobulk_preprocessing_SE = pseudobulk_preprocessing(reference_label_fine,
                                                       preprocessing_output_S, 
                                                       !!sample_column)

test_that("input_read_RNA_assay_works", {
  expect_s4_class(input_read_RNA_assay, "Seurat")
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
  expect_true(nrow(empty_droplets_tbl) < nrow(input_read_RNA_assay))
})

test_that("cell_cycle_score_works", {
  expect_s3_class(cell_cycle_score_tbl, "tbl_df")
  expected_colnames <- c("S.Score", "G2M.Score", "Phase")
  expect_true(all(expected_colnames %in% colnames(cell_cycle_score_tbl)))
})

test_that("annotation_label_transfer_works", {
  
  if (!is.null(reference_azimuth)) {
    # Expect the output to be a tibble
    expect_equal(ncol(annotation_label_transfer_tbl), 10)
  } else {
    expect_equal(ncol(annotation_label_transfer_tbl), 9)
  }
})

test_that("alive_identification_works", {
  expect_s3_class(alive_identification_tbl, "tbl_df")
  expected_colnames <- c("subsets_Mito_sum", "subsets_Mito_detected", "subsets_Mito_percent")
  expect_true(all(expected_colnames %in% colnames(alive_identification_tbl)))
})

test_that("non_batch_variation_removal_S_dimensions", {
  num_features_input = nrow(input_read_RNA_assay@assays$RNA@counts)
  num_cells_input = ncol(input_read_RNA_assay@assays$RNA@counts)

  num_features_non_batch = nrow(non_batch_variation_removal_S@assays$SCT@counts)
  num_cells_non_batch = ncol(non_batch_variation_removal_S@assays$SCT@counts)
  # Expect less features 
  expect_true(num_features_non_batch < num_features_input)
  # Expect less cells 
  expect_true(num_cells_non_batch < num_cells_input)
})

test_that("Doublet_identification_works", {
  # Expect a tibble
  expect_s3_class(doublet_identification_tbl, "tbl_df")
  
  expected_colnames <- c("scDblFinder.class", "scDblFinder.score", "scDblFinder.weighted", "scDblFinder.cxds_score")
  expect_true(all(expected_colnames %in% colnames(doublet_identification_tbl)))
})


test_that("Preprocessing_works", {
  expect_s4_class(preprocessing_output_S, "Seurat")
})




