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
})

test_that("cell_cycle_score_works", {
  expect_s3_class(cell_cycle_score_tbl, "tbl_df")
})

test_that("annotation_label_transfer_works", {
  
  # If seurat is provided
  expect_s3_class(annotation_label_transfer_tbl, "tbl_df")
  # If seurat is not
  expect_s3_class(annotation_label_transfer_tbl, "tbl_df")
})

test_that("alive_identification_works", {
  expect_s3_class(alive_identification_tbl, "tbl_df")
})

test_that("non-batch_variation_removal_works", {
  expect_s4_class(non_batch_variation_removal_S, "Seurat")
})

test_that("Doublet_identification_works", {
  expect_s3_class(doublet_identification_tbl, "tbl_df")
})

test_that("Preprocessing_works", {
  expect_s4_class(preprocessing_output_S, "Seurat")
})




