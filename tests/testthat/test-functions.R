## Load packages
# my_packages <- c( "HPCell",
#                   "glue",
#                   "readr",
#                   "dplyr",
#                   "tidyr",
#                   "ggplot2",
#                   "purrr",
#                   "Seurat",
#                   "tidyseurat",
#                   "glue",
#                   "scater",
#                   "DropletUtils",
#                   "EnsDb.Hsapiens.v86",
#                   "here",
#                   "stringr",
#                   "readr",
#                   "rlang",
#                   "scuttle",
#                   "scDblFinder",
#                   "ggupset",
#                   "tidySummarizedExperiment",
#                   "broom",
#                   "tarchetypes",
#                   "SeuratObject",
#                   "SingleCellExperiment",
#                   "SingleR",
#                   "celldex",
#                   "tidySingleCellExperiment",
#                   "tibble",
#                   "magrittr",
#                   "tidybulk", 
#                   "scRNAseq")
# lapply(my_packages, require, character.only = TRUE)

### Load input data 
file_path = "input_data.rds"
if(!file.exists(file_path)) 
{
  single_cell_data = scRNAseq::HeOrganAtlasData(ensembl=FALSE,location=FALSE)[,1:100]
  single_cell_data |> Seurat::as.Seurat(data = NULL) |> saveRDS(file_path)
}
###Load reference data 
input_reference_path <- "reference_azimuth.rds"
if(!file.exists(input_reference_path)) 
  {
  reference_url<- "https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat"
  download.file(reference_url, input_reference_path, cacheOK = TRUE)
  LoadH5Seurat(input_reference_path) |> saveRDS(input_reference_path)
}

## Define arguments 
filtered <- "TRUE"
tissue <- "pbmc"
RNA_assay_name<- "originalexp"
input_file<- readRDS(file_path)
reference_azimuth = readRDS(input_reference_path)

# Function testing
# test_that("add_RNA_assay_works", {
#   RNA_assay = add_RNA_assay(readRDS(file_path), RNA_assay_name)
#   expect_s4_class(res, "Seurat")})
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



test_that("empty_droplets_works", {
  expect_s3_class(empty_droplets_tbl, "tbl_df")
})

test_that("cell_cycle_score_works", {
  expect_s3_class(cell_cycle_scoring, "tbl_df")
})

test_that("annotation_label_transfer_works", {
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



# test_that("pseudobulk_preprocessing_works", {
#   expect_s4_class(res$pseudobulk_by_sample, "SummarizedExperiment")
#   expect_s4_class(res$pseudobulk_by_sample_and_cell_type, "SummarizedExperiment")
# })

# test_that("reference_label_id_works", {
#   res = reference_label_id(tissue, reference_label)
#   expect_type(res, "character")
# })
# 

# 
# 
# se = 
#   tidySummarizedExperiment::se |> 
#   tidybulk::keep_abundant()
# 
# test_that("de analysis single", {
#   
#   se |> 
#     hpcell_test_differential_abundance(~ dex + (1 | cell))
# })
# 
# test_that("de analysis multi", {
#   
# 
#   tibble(name = c("my_data_1", "my_data_2"), data = list(se, se)) |> 
#       hpcell_map_test_differential_abundance(~ dex + (1 | cell), data)
#   
# 
# })
# 
# 
# test_that("split", {
#   
#   tibble(name = c("my_data_1"), data = list(se)) |> 
#     mutate(split = 5) |> 
#     map_split_se_by_gene(data, split) |> 
#     nrow() |> 
#     expect_equal(5)
# 
# })



# alive_identification_B14<- readRDS("/home/users/allstaff/si.j/test_HPCell/results/preprocessing_results/alive_identification/CB150T04X__batch14__alive_identification_output.rds")
# 
# doublet_identification_B14<- readRDS("/home/users/allstaff/si.j/test_HPCell/results/preprocessing_results/doublet_identification/CB150T04X__batch14__doublet_identification_output.rds")
# 
# annotation_label_transfer_B14<- readRDS("/home/users/allstaff/si.j/test_HPCell/results/preprocessing_results/annotation_label_transfer/CB150T04X__batch14__annotation_label_transfer_output.rds")
# 
# non_batch_variation_removal <- readRDS("/home/users/allstaff/si.j/test_HPCell/results/preprocessing_results/non_batch_variation_removal/CB150T04X__batch14__non_batch_variation_removal_output.rds")
# 
# cell_cycle_scoring <- readRDS("/home/users/allstaff/si.j/test_HPCell/results/preprocessing_results/cell_cycle_scoring/CB150T04X__batch14__cell_cycle_scoring_output.rds")
# 
# reference_azimuth <- readRDS("/home/users/allstaff/si.j/HPCell_ARCHIVE/data/Data/jiayi_files/reference_azimuth.rds")


