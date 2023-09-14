# Loading packages
my_packages <- c( "jascap",
                  "glue",
                  "readr",
                  "dplyr",
                  "tidyr",
                  "ggplot2",
                  "purrr",
                  "Seurat",
                  "tidyseurat",
                  "glue",
                  "scater",
                  "DropletUtils",
                  "EnsDb.Hsapiens.v86",
                  "here",
                  "stringr",
                  "readr",
                  "rlang",
                  "scuttle",
                  "scDblFinder",
                  "ggupset",
                  "tidySummarizedExperiment",
                  "broom",
                  "tarchetypes",
                  "SeuratObject",
                  "SingleCellExperiment",
                  "SingleR",
                  "celldex",
                  "tidySingleCellExperiment",
                  "tibble",
                  "magrittr",
                  "tidybulk")
lapply(my_packages, require, character.only = TRUE)
reference_label<-"none"
reference_label_fine<- "monaco_first.labels.fine"
input_path_preprocessing_output <- list(readRDS("~/test_jascap/results/preprocessing_results/preprocessing_output/CB150T04X__batch14__preprocessing_output_output.rds"), 
                                        readRDS("~/test_jascap/results/preprocessing_results/preprocessing_output/CB291T01X__batch8__preprocessing_output_output.rds"))
#
# #Defining inputs
input_file_B14<- readRDS('/home/users/allstaff/si.j/test_jascap/input/CB150T04X__batch14.rds')

empty_droplets_tbl_B14<- readRDS("/home/users/allstaff/si.j/test_jascap/results/preprocessing_results/empty_droplet_identification/CB150T04X__batch14__empty_droplet_identification_output.rds")

alive_identification_B14<- readRDS("/home/users/allstaff/si.j/test_jascap/results/preprocessing_results/alive_identification/CB150T04X__batch14__alive_identification_output.rds")

doublet_identification_B14<- readRDS("/home/users/allstaff/si.j/test_jascap/results/preprocessing_results/doublet_identification/CB150T04X__batch14__doublet_identification_output.rds")

annotation_label_transfer_B14<- readRDS("/home/users/allstaff/si.j/test_jascap/results/preprocessing_results/annotation_label_transfer/CB150T04X__batch14__annotation_label_transfer_output.rds")

non_batch_variation_removal <- readRDS("/home/users/allstaff/si.j/test_jascap/results/preprocessing_results/non_batch_variation_removal/CB150T04X__batch14__non_batch_variation_removal_output.rds")

cell_cycle_scoring <- readRDS("/home/users/allstaff/si.j/test_jascap/results/preprocessing_results/cell_cycle_scoring/CB150T04X__batch14__cell_cycle_scoring_output.rds")

reference_azimuth <- readRDS("/home/users/allstaff/si.j/jascap_ARCHIVE/data/Data/jiayi_files/reference_azimuth.rds")
filtered <- "filtered"
tissue <- "pbmc"
# Function testing
test_that("empty_droplets_works", {
  res = empty_droplet_id(input_file_B14,
                         filtered)
  expect_s3_class(res, "tbl_df")
})

# test_that("cell_cycle_score_works", {
#   res = cell_cycle_scoring(input_file_B14,
#                            empty_droplets_tbl_B14)
#   expect_s3_class(res, "tbl_df")
# })
# 
test_that("annotation_label_transfer_works", {
  res = annotation_label_transfer(input_file_B14,
                                  reference_azimuth,
                                  empty_droplets_tbl_B14)
  expect_s3_class(res, "tbl_df")
})

# test_that("alive_identification_works", {
#   res = alive_identification(input_file_B14,
#                              empty_droplets_tbl_B14,
#                              annotation_label_transfer_B14
#   )
#   expect_s3_class(res, "tbl_df")
# })
# 
# test_that("non-batch_variation_removal_works", {
#   res = non_batch_variation_removal(input_file_B14,
#                                     empty_droplets_tbl_B14,
#                                     alive_identification_B14,
#                                     cell_cycle_scoring)
#   expect_s4_class(res, "Seurat")
#   
# })
# 
# test_that("Doublet_identification_works", {
#   res = doublet_identification(input_file_B14,
#                                empty_droplets_tbl_B14,
#                                alive_identification_B14,
#                                annotation_label_transfer_B14,
#                                reference_label_fine)
#   expect_s3_class(res, "tbl_df")
# })
# test_that("Preprocessing_works", {
#   res = preprocessing_output(tissue,
#                              non_batch_variation_removal,
#                              alive_identification_B14,
#                              cell_cycle_scoring,
#                              annotation_label_transfer_B14,
#                              doublet_identification_B14)
#   expect_s4_class(res, "Seurat")
# })
# 
# test_that("pseudobulk_preprocessing_works", {
#   res = pseudobulk_preprocessing(
#     reference_label_fine, input_path_preprocessing_output
#   )
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
