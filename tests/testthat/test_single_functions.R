library(testthat)
library(HPCell)
library(Seurat)
library(scRNAseq)
## Define arguments 
filter_empty_droplets <- "TRUE"
tissue <- "pbmc"
RNA_assay_name<- "originalexp"
#reference_azimuth<- NULL

input_seurat_abc = 
  HeOrganAtlasData(ensembl=FALSE,location=FALSE)|> 
  as.Seurat(data = NULL) 

sample_column<- "Tissue"
## Defining functions 

reference_label_fine = HPCell:::reference_label_fine_id(tissue)
empty_droplets_tbl = HPCell:::empty_droplet_id(input_seurat_abc, filter_empty_droplets = TRUE)

# Define output from annotation_label_transfer 
annotation_label_transfer_tbl = HPCell:::annotation_label_transfer(input_seurat_abc,
                                                          empty_droplets_tbl)

# Define output from alive_identification
alive_identification_tbl = HPCell:::alive_identification(input_seurat_abc,
                                                empty_droplets_tbl,
                                                annotation_label_transfer_tbl)


# Define output from doublet_identification
doublet_identification_tbl = HPCell:::doublet_identification(input_seurat_abc,
                                                    empty_droplets_tbl,
                                                    alive_identification_tbl,
                                                    annotation_label_transfer_tbl,
                                                    reference_label_fine)

# Define output from cell_cycle_scoring
cell_cycle_score_tbl = HPCell:::cell_cycle_scoring(input_seurat_abc, empty_droplets_tbl)

# Define output from non_batch_variation_removal
non_batch_variation_removal_S = non_batch_variation_removal(input_seurat_abc,
                                                            empty_droplets_tbl,
                                                            alive_identification_tbl,
                                                            cell_cycle_score_tbl)
# Define output from preprocessing_output
preprocessing_output_S = HPCell:::preprocessing_output(tissue,
                                              non_batch_variation_removal_S,
                                              alive_identification_tbl,
                                              cell_cycle_score_tbl,
                                              annotation_label_transfer_tbl,
                                              doublet_identification_tbl)

# Define output from pseudobulk_preprocessing
pseudobulk_preprocessing_SE = HPCell:::pseudobulk_preprocessing(reference_label_fine, 
                                                       preprocessing_output_S, 
                                                       sample_column)

create_pseudobulk_sample = HPCell:::create_pseudobulk(preprocessing_output_S, assays = "RNA", x = c(Tissue, Cell_type_in_each_tissue))

pseudobulk_merge_all_samples = pseudobulk_merge(list(create_pseudobulk_sample), assays = "RNA", x = c(Tissue))

# Testing function outputs are as expected 
test_that("input_seurat_works", {
  expect_s4_class(input_seurat, "Seurat")
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
  expect_true(nrow(empty_droplets_tbl) < nrow(input_seurat))
})

test_that("cell_cycle_score_works", {
  expect_s3_class(cell_cycle_score_tbl, "tbl_df")
  expected_colnames <- c("S.Score", "G2M.Score", "Phase")
  expect_true(all(expected_colnames %in% colnames(cell_cycle_score_tbl)))
})

test_that("annotation_label_transfer_works", {
  
  # if (!is.null(reference_azimuth)) {
  #   # Expect the output to be a tibble
  #   expect_equal(ncol(annotation_label_transfer_tbl), 10)
  # } else {
    expect_equal(ncol(annotation_label_transfer_tbl), 9)
  # }
})

test_that("alive_identification_works", {
  expect_s3_class(alive_identification_tbl, "tbl_df")
  expected_colnames <- c("subsets_Mito_sum", "subsets_Mito_detected", "subsets_Mito_percent")
  expect_true(all(expected_colnames %in% colnames(alive_identification_tbl)))
})

test_that("non_batch_variation_removal_S_dimensions", {
  num_features_input = nrow(input_seurat)
  num_cells_input = ncol(input_seurat)

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

# test_that("pseudobulk_preprocessing handles input lists", {
#   expect_s4_class(pseudobulk_preprocessing_SE[[1]], "SummarizedExperiment") 
#   expect_s4_class(pseudobulk_preprocessing_SE[[2]], "SummarizedExperiment")                
# })

test_that("pseudobulk_sample_works", {
  expect_
})

test_that("pseudobulk_preprocessing_works", {
  expect_s4_class(pseudobulk_preprocessing_SE, "Seurat")
})

## Test reports 
# Subset 2 tissues to test reports
subset_tissue <- subset(input_seurat_abc, subset = Tissue %in% c("Heart", "Trachea"))
input_seurat<- add_RNA_assay(subset_tissue, "originalexp")
heart <- subset(input_seurat, subset = Tissue == "Heart")
trachea <- subset(input_seurat, subset = Tissue == "Trachea")
input_seurat_list <- c(heart, trachea)


unique_idents <- unique(input_seurat$orig.ident)
list_of_subsets <- lapply(unique_idents, function(ident) {
  subset(input_seurat, subset = orig.ident == ident)
})

# Test empty droplets
empty_droplets_tissue_list <- lapply(input_seurat_list, function(df) {
  HPCell:::empty_droplet_id(df, filter_empty_droplets = TRUE)
})

annotation_label_transfer_tbl_list <- mapply(FUN = HPCell:::annotation_label_transfer, 
                                             input_seurat_list, 
                                             empty_droplets_tissue_list,
                                             SIMPLIFY = FALSE)

alive_identification_tbl_list <- mapply(FUN = HPCell:::alive_identification, 
                                        input_seurat_list, 
                                        empty_droplets_tissue_list, annotation_label_transfer_tbl_list,
                                        SIMPLIFY = FALSE)

doublet_identification_tbl_list <- mapply(FUN = HPCell:::doublet_identification, 
                                          input_seurat_list, 
                                          empty_droplets_tissue_list, 
                                          alive_identification_tbl_list, 
                                          annotation_label_transfer_tbl_list, 
                                          reference_label_fine
)

# Unit test 
test_that("R Markdown render works", {
  # Define output paths
  input_path <- paste0(system.file(package = "HPCell"), "/rmd/Empty_droplet_report.Rmd")
  output_path <- paste0(system.file(package = "HPCell"), "/Empty_droplet_report.html")
  
  # Test execution: Render the R Markdown file
  rmarkdown::render(
    input = input_path,
    output_file = output_path,
    params = list(x1 = input_seurat_list, x2 = empty_droplets_tissue_list)
  )
  
  # Assertions
  expect_true(file.exists(output_path), info = "Output file should exist")
  # Add more assertions as needed, e.g., checking file content, format, etc.
})

preprocessing_output_S




