library(testthat)
library(HPCell)
library(Seurat)
library(scRNAseq)
## Define arguments 
filter_empty_droplets <- "TRUE"
tissue <- "pbmc"
RNA_assay_name<- "originalexp"

input_seurat_abc = 
  HeOrganAtlasData(ensembl=FALSE,location=FALSE)|> 
  as.Seurat(data = NULL) |>
  subset(subset = Tissue %in% c("Blood")) 

cell_type_column <- "Cell_type_in_each_tissue"
# sample_column<- "Tissue"
## Defining functions 
# 
reference_label_fine = HPCell:::reference_label_fine_id(tissue)

empty_droplets_tbl = HPCell:::empty_droplet_threshold(input_seurat_abc, 
                                                      feature_nomenclature = "symbol")

# Define output from annotation_label_transfer 
annotation_label_transfer_tbl = HPCell:::annotation_label_transfer(input_seurat_abc,
                                                          empty_droplets_tbl,
                                                          reference_azimuth = "pbmcref",
                                                          feature_nomenclature = "symbol")

# Define output from cell_type_ensembl_harmonised
cell_type_ensemble_tbl = HPCell:::cell_type_ensembl_harmonised(input_seurat_abc, 
                                                              annotation_label_transfer_tbl)

# Define output from alive_identification
alive_identification_tbl = HPCell:::alive_identification(input_seurat_abc,
                                                empty_droplets_tbl,
                                                cell_type_ensemble_tbl,
                                                cell_type_column = "cell_type_unified_ensemble",
                                                feature_nomenclature = "symbol")


# Define output from doublet_identification
doublet_identification_tbl = HPCell:::doublet_identification(input_seurat_abc,
                                                    empty_droplets_tbl,
                                                    alive_identification_tbl,
                                                    cell_type_ensemble_tbl,
                                                    reference_label_fine = "cell_type_unified_ensemble")

# Define output from cell_cycle_scoring
cell_cycle_score_tbl = HPCell:::cell_cycle_scoring(input_seurat_abc, empty_droplets_tbl)

# Define output from non_batch_variation_removal
non_batch_variation_removal_S = HPCell:::non_batch_variation_removal(input_seurat_abc,
                                                            empty_droplets_tbl,
                                                            alive_identification_tbl,
                                                            cell_cycle_score_tbl, 
                                                            assay = NULL)
# Define output from preprocessing_output
preprocessing_output_S = HPCell:::preprocessing_output(tissue,
                                              non_batch_variation_removal_S,
                                              alive_identification_tbl,
                                              cell_cycle_score_tbl,
                                              annotation_label_transfer_tbl,
                                              doublet_identification_tbl)

# Calculate metacell for a sample cell type
metacell_per_cell_type <- HPCell:::calculate_metacell_for_a_sample_per_cell_type(input_seurat_abc,
                                                                                 min_cells_per_metacell = 10)

# Calculate metacell membership
metacell_tbl <- split_sample_cell_type_calculate_metacell_membership(input_seurat_abc, 
                                                                     input_seurat_abc[[]] |> 
                                                                       rownames_to_column(var = ".cell") |> 
                                                                       as_tibble(),
                                                                     empty_droplets_tbl,
                                                                     alive_identification_tbl,
                                                                     doublet_identification_tbl,
                                                                     x = cell_type_column,
                                                                     min_cells_per_metacell = 10)

# Define output from cell_communication
cell_communication_tbl = HPCell:::cell_communication(input_seurat_abc,
                                                     empty_droplets_tbl = NULL,
                                                     alive_identification_tbl = NULL,
                                                     doublet_identification_tbl = NULL,
                                                     cell_type_tbl = input_seurat_abc[[]] |> 
                                                       rownames_to_column(var = ".cell") |> 
                                                       as_tibble() |> mutate(sample_id = "sample1"),
                                                     assay = NULL,
                                                     cell_type_column = "Cell_type_in_each_tissue",
                                                     feature_nomenclature = "symbol")

# empty_droplets_tbl = HPCell:::empty_droplet_id(input_seurat_list[[1]], filter_empty_droplets = TRUE)
# 
# # Define output from annotation_label_transfer 
# annotation_label_transfer_tbl = HPCell:::annotation_label_transfer(input_seurat_list[[1]],
#                                                           empty_droplets_tbl)
# 
# # Define output from alive_identification
# alive_identification_tbl = HPCell:::alive_identification(input_seurat_list[[1]],
#                                                 empty_droplets_tbl,
#                                                 annotation_label_transfer_tbl, 
#                                                 assay = NULL)
# 
# 
# # Define output from doublet_identification
# doublet_identification_tbl = HPCell:::doublet_identification(input_seurat_list[[1]],
#                                                     empty_droplets_tbl,
#                                                     alive_identification_tbl,
#                                                     annotation_label_transfer_tbl,
#                                                     reference_label_fine)
# 
# # Define output from cell_cycle_scoring
# cell_cycle_score_tbl = HPCell:::cell_cycle_scoring(input_seurat_list[[1]], empty_droplets_tbl)
# 
# # Define output from non_batch_variation_removal
# non_batch_variation_removal_S = HPCell:::non_batch_variation_removal(input_seurat_list[[1]],
#                                                             empty_droplets_tbl,
#                                                             alive_identification_tbl,
#                                                             cell_cycle_score_tbl, 
#                                                             assay = NULL)
# # Define output from preprocessing_output
# preprocessing_output_S = HPCell:::preprocessing_output(tissue,
#                                               non_batch_variation_removal_S,
#                                               alive_identification_tbl,
#                                               cell_cycle_score_tbl,
#                                               annotation_label_transfer_tbl,
#                                               doublet_identification_tbl)

# Define output from pseudobulk_preprocessing
# pseudobulk_preprocessing_SE = HPCell:::pseudobulk_preprocessing(reference_label_fine, 
#                                                        preprocessing_output_S, 
#                                                        sample_column)
# 
# create_pseudobulk_sample = HPCell::create_pseudobulk(preprocessing_output_S, assays = "SCT", x = c(Tissue, Cell_type_in_each_tissue))

# For a list of preprocessing outputs
# create_pseudobulk_sample_list = mapply(FUN = create_pseudobulk, 
#                                        preprocessing_output_S_list, 
#                                        assays = "RNA", 
#                                        x = c(Tissue, Cell_type_in_each_tissue))

create_pseudobulk_sample_list <- lapply(preprocessing_output_S_list, function(obj) {
  create_pseudobulk(obj, assays = "originalexp", x = c(Tissue, Cell_type_in_each_tissue))
})

pseudobulk_merge_all_samples = pseudobulk_merge(create_pseudobulk_sample_list, assays = "RNA", x = c(Tissue))

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


unique_idents <- 
  c(
    get_unique_tissues(input_seurat_list[[1]]),
    get_unique_tissues(input_seurat_list[[2]])
  )


if(is.null(assay)) assay = input_seurat_abc@assays |> names() |> extract2(1)
#reference_azimuth<- NULL
process_seurat_object <- function(input_seurat_abc, assay = NULL) {
  if(is.null(assay)) {
    assay <- input_seurat_abc@assays |> names() |> extract2(1)
  }
  # Return the updated assay (or the original one if it was not NULL)
  return(assay)
}

assay<- process_seurat_object(input_seurat_list[[1]])

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
                                          reference_label_fine, 
                                          SIMPLIFY = FALSE)

cell_cycle_score_list <- mapply(FUN = HPCell:::cell_cycle_scoring, 
                                input_seurat_list, 
                                empty_droplets_tissue_list,
                                SIMPLIFY = FALSE)

non_batch_variation_removal_list = mapply(FUN = HPCell:::non_batch_variation_removal, 
                                          input_seurat_list,
                                          empty_droplets_tissue_list,
                                          alive_identification_tbl_list,
                                          cell_cycle_score_list,
                                          assay = assay,
                                          SIMPLIFY = FALSE)

preprocessing_output_S_list = mapply(FUN = HPCell:::preprocessing_output, 
                                     tissue,
                                     non_batch_variation_removal_list,
                                     alive_identification_tbl_list,
                                     cell_cycle_score_list,
                                     annotation_label_transfer_tbl_list,
                                     doublet_identification_tbl_list)


create_pseudobulk_sample_list = mapply(FUN = create_pseudobulk, 
                                       preprocessing_output_S_list, 
                                       assays = assay, 
                                       x = c(Tissue, Cell_type_in_each_tissue))

create_pseudobulk_sample_heart<- create_pseudobulk(preprocessing_output_S_list[[1]], assays = assay, x = c(Tissue, Cell_type_in_each_tissue))
create_pseudobulk_sample_trachea <- create_pseudobulk(preprocessing_output_S_list[[2]], assays = assay, x = c(Tissue, Cell_type_in_each_tissue))

## Fibrosis 
create_pseudobulk_sample_list <- lapply(preprocessing_output_S_list, function(obj) {
  create_pseudobulk(obj, assays = "RNA", x = c(sampleName, cellAnno))
})
pseudobulk_merge_all_samples = pseudobulk_merge(create_pseudobulk_sample_list, assays = "RNA", x = c(Tissue))
# create_pseudobulk_sample_list = mapply(FUN = create_pseudobulk,
#                                        preprocessing_output_S_list,
#                                        assays = assay,
#                                        x = c(sampleName, Cell_type_in_each_tissue))
 
# create_pseudobulk_sample_heart<- create_pseudobulk(preprocessing_output_S_list[[1]], assays = assay, x = c(Tissue, Cell_type_in_each_tissue))
# create_pseudobulk_sample_trachea <- create_pseudobulk(preprocessing_output_S_list[[2]], assays = assay, x = c(Tissue, Cell_type_in_each_tissue))
# 
# create_pseudobulk_sample_list<- list(create_pseudobulk_sample_heart, create_pseudobulk_sample_trachea)

# create_pseudobulk_sample_list <- lapply(preprocessing_output_S_list, function(obj) {
#   create_pseudobulk(obj, assays = NULL, x = c(Tissue, Cell_type_in_each_tissue))
# })
# 
# pseudobulk_merged_results <- pseudobulk_merge(create_pseudobulk_sample_list, assays, x)

# Test calc_UMAP_dbl_report
# calc_UMAP_result<- calc_UMAP(input_seurat)

calc_UMAP_result_list<- lapply(input_seurat_list, function(df) {
  #HPCell:::calc_UMAP(df)
  calc_UMAP(df, assay = NULL)
})

# Unit test 
test_that("R Markdown render empty droplet works", {
  # Define output paths
  input_path <- paste0(system.file(package = "HPCell"), "/rmd/Empty_droplet_report.Rmd")
  output_path <- paste0(system.file(package = "HPCell"), "/Empty_droplet_report.html")
  
  # Test execution: Render the R Markdown file
  rmarkdown::render(
    input = input_path,
    output_file = output_path,
    params = list(x1 = input_seurat_list, 
                  x2 = empty_droplets_tissue_list, 
                  x3 = annotation_label_transfer_tbl_list, 
                  x4 = unique_idents, 
                  x5 = sample_column)
  )
  
  # Assertions
  expect_true(file.exists(output_path), info = "Output file should exist")
  # Add more assertions as needed, e.g., checking file content, format, etc.
})

test_that("calc_UMAP returns correctly structured tibble", {
  # Assuming calc_UMAP is your function and input_seurat is the input
  output <- calc_UMAP(input_seurat)
  
  # Check if output is a tibble
  expect_true(is_tibble(output))
})

test_that("R Markdown render doublet identification works", {
  input_path <- paste0(system.file(package = "HPCell"), "/rmd/Doublet_identification_report.Rmd")
  output_path <- paste0(system.file(package = "HPCell"), "/Doublet_identification_report.html")
  
  rmarkdown::render(
    input = input_path,
    output_file = output_path,
    params = list(x1 = input_seurat_list,
                  x2 = calc_UMAP_result_list,
                  x3 = doublet_identification_tbl_list,
                  x4 = annotation_label_transfer_tbl_list, 
                  x5 = sample_column |> enquo(), 
                  x6 = cell_type_annotation_column |> enquo())
    )
  expect_true(file.exists(output_path), info = "Output file should exist")
})

test_that("R Markdown render pseudobulk analysis works", {
  input_path <- paste0(system.file(package = "HPCell"), "/rmd/pseudobulk_analysis_report.Rmd")
  output_path <- paste0(system.file(package = "HPCell"), "/Pseudobulk_analysis_report.html")
  rmarkdown::render(
    input = input_path,
    output_file = output_path,
    params = list(x1 = pseudobulk_merge_all_samples)
  )
})


## Technical_variation_report 

rmarkdown::render(
  input = paste0(system.file(package = "HPCell"), "/rmd/Technical_variation_report.Rmd"),
  output_file = paste0(system.file(package = "HPCell"), "/Technical_variation_report.html"),
  params = list(
    x1 = input_seurat_list,
    x2 = empty_droplets_tissue_list)
)

path<- paste0(system.file(package = "HPCell"), "extdata/Test.Rmd")

## Testing in Targets 

file_paths <- tar_read(read_file_files, store = store)
data_container_type <- tar_read(data_container_type_file, store = store)

input_file <- list()

# Loop over elements in file_paths list
for (i in seq_along(file_paths)) {
  input_file[[i]] <- read_data_container(file_paths[[i]], container_type = data_container_type)
}

## Empty Droplets
rmarkdown::render(
  input =  paste0(system.file(package = "HPCell"), "/rmd/Empty_droplet_report.Rmd"),
  output_file = paste0(system.file(package = "HPCell"), "/Empty_droplet_report.html"),
  params = list(x1 = read_data_container(tar_read(file_path, store = store), container_type = tar_read(data_container_type_file, store = store)), 
                x2 = tar_read(empty_droplets_tbl, store = store),
                x3 = tar_read(annotation_label_transfer_tbl, store = store),
                x4 = tar_read(unique_tissues, store = store), 
                x5 = tar_read(sample_column, store = store)|> quo_name())
)

## Doublet identification 
rmarkdown::render(
  input = paste0(system.file(package = "HPCell"), "/rmd/Doublet_identification_report.Rmd"),
  output_file = paste0(system.file(package = "HPCell"),"/Doublet_identification_report.html"),
  params = list(x1 = tar_read(input_read, store = store),
                x2 = tar_read(calc_UMAP_dbl_report, store = store),
                x3 = tar_read(doublet_identification_tbl, store = store),
                x4 = tar_read(annotation_label_transfer_tbl, store = store), 
                x5 = tar_read(sample_column, store = store) |> quo_name(), 
                x6 = tar_read(cell_type_annotation_column, store = store) |> quo_name()
))

## Technical variation 
rmarkdown::render(
  input = paste0(system.file(package = "HPCell"), "/rmd/Technical_variation_report.Rmd"),
  output_file = paste0(system.file(package = "HPCell"), "/Technical_variation_report.html"),
  params = list(x1 = tar_read(input_read, store = store),
                x2 = tar_read(empty_droplets_tbl, store = store), 
                x3 = tar_read(variable_gene_list, store = store), 
                x4 = tar_read(calc_UMAP_dbl_report, store = store), 
                x5 = tar_read(sample_column, store = store) |> quo_name()
  )
)

## Pseudobulk analysis report 
rmarkdown::render(
  input = paste0(system.file(package = "HPCell"), "/rmd/pseudobulk_analysis_report.Rmd"),
  output_file = paste0(system.file(package = "HPCell"), "/pseudobulk_analysis_report.html"),
  params = list(x1 = tar_read(pseudobulk_merge_all_samples, store = store), 
                x2 = tar_read(sample_column, store = store) |> quo_name(), 
                x3 = tar_read(cell_type_annotation_column, store = store) |> quo_name())
)


# 
# tissues <- unique(input_seurat$Tissue)
# seurat_objects <- list()
# 
# for (tissue in tissues) {
#   seurat_objects[[tissue]] <- subset(input_seurat, idents = tissue)
# }
# seurat_objects[["Heart"]]
# seurat_objects[["Trachea"]]


## Testing fibrosis dataset 
input_a<- readRDS("~/HPCell/fibrosis_data/GSE122960___GSM3489182.rds")
input_b <- readRDS("~/HPCell/fibrosis_data/GSE135893_cHP___THD0001.rds")
input_seurat_list<- c(input_a, input_b)
sample_column<- "sampleName"
cell_type_annotation_column <- "cellAnno"
filter_empty_droplets <- "TRUE"
tissue <- "pbmc"
reference_label_fine = HPCell:::reference_label_fine_id(tissue)

if(is.null(assay)) assay = input_a@assays |> names() |> extract2(1)
#reference_azimuth<- NULL
process_seurat_object <- function(input_a, assay = NULL) {
  if(is.null(assay)) {
    assay <- input_a@assays |> names() |> extract2(1)
  }
  # Return the updated assay (or the original one if it was not NULL)
  return(assay)
}

assay<- process_seurat_object(input_seurat_list[[1]])

# library(SeuratData)
# # devtools::install_github("satijalab/azimuth")
# library(Azimuth)
# library(Seurat)
# options(Seurat.object.assay.version = "v5")
# obj <- LoadData("pbmcsca") |>
#   UpdateSeuratObject() |>
#   RunAzimuth(reference = "pbmcref") |>
#   SCTransform()
# obj |> saveRDS("dev/reference_azimuth.rds")


library(HPCell)
library(crew)
library(crew.cluster)

# library(magrittr)
# library(tidySingleCellExperiment)
# library(Seurat)
# library(SeuratData)
# InstallData("pbmc3k")
# options(Seurat.object.assay.version = "v5")
input_seurat <-
  LoadData("pbmc3k") |>
  _[,1:500]
# 
# change_seurat_counts = function(data){
# 
#   data@assays$RNA$counts = data@assays$RNA$counts * runif(min = 0, max = 2, n = length(data@assays$RNA$counts))
#   data@assays$RNA$data = data@assays$RNA$counts
#   data
# }
# input_seurat |> mutate(condition = "treated") |> change_seurat_counts() |> mutate(group = 1) |>   saveRDS("dev/input_seurat_treated_1.rds")
# input_seurat |> mutate(condition = "treated") |> change_seurat_counts() |>  mutate(group = 2) |> saveRDS("dev/input_seurat_treated_2.rds")
# input_seurat |> mutate(condition = "untreated") |> change_seurat_counts() |>  mutate(group = 1) |> saveRDS("dev/input_seurat_UNtreated_1.rds")
# input_seurat |> mutate(condition = "untreated") |> change_seurat_counts() |>  mutate(group = 2) |> saveRDS("dev/input_seurat_UNtreated_2.rds")

# input_seurat |> mutate(condition = "treated") |> change_seurat_counts() |> as.SingleCellExperiment() |>   saveRDS("dev/input_seurat_treated_1_SCE.rds")
# input_seurat |> mutate(condition = "treated") |> change_seurat_counts() |> as.SingleCellExperiment() |>  saveRDS("dev/input_seurat_treated_2_SCE.rds")
# input_seurat |> mutate(condition = "untreated") |> change_seurat_counts() |> as.SingleCellExperiment() |>  saveRDS("dev/input_seurat_UNtreated_1_SCE.rds")
# input_seurat |> mutate(condition = "untreated") |> change_seurat_counts() |> as.SingleCellExperiment() |>  saveRDS("dev/input_seurat_UNtreated_2_SCE.rds")

# library(SeuratData)
# InstallData("pbmcsca")
# pbmcsca <- LoadData("pbmcsca") # save this to disk, so you can recall every time you execute HPCell

computing_resources = crew_controller_local(workers = 8) #resource_tuned_slurm

# tier = rep(c("tier_1", "tier_2"), times = 6),
# computing_resources = list(
#   
#   crew_controller_local(
#     name = "tier_1",
#     workers = 4
#   ),
#   crew_controller_local(
#     name = "tier_2",
#     workers = 4
#   )
# )

  computing_resources = list(

  crew_controller_slurm(
    name = "tier_1",
    slurm_memory_gigabytes_per_cpu = 5,
    slurm_cpus_per_task = 1,
    workers = 50,
    tasks_max = 5,
    verbose = T, 
    seconds_idle = 30
  ),
  crew_controller_slurm(
    name = "tier_2",
    slurm_memory_gigabytes_per_cpu = 10,
    slurm_cpus_per_task = 1,
    workers = 50,
    tasks_max = 5,
    verbose = T, 
    seconds_idle = 30
  )
)

#  Slurm resources
# computing_resources =
#   crew.cluster::crew_controller_slurm(
#     slurm_memory_gigabytes_per_cpu = 5,
#     workers = 500,
#     tasks_max = 5,
#     verbose = T,
#     slurm_cpus_per_task = 1
#   )

{ foo_function = function(x, y){x |> dplyr::mutate(foo = y)} } |> 
  substitute() |> 
  deparse() |> 
  readr::write_lines("dev/my_custom_script.R")

# # Define and execute the pipeline
file_list = 
#   c("dev/input_seurat_treated_1.rds",
#   "dev/input_seurat_treated_2.rds",
#   "dev/input_seurat_UNtreated_1.rds",
#   "dev/input_seurat_UNtreated_2.rds") |>
#   purrr::map_chr(here::here) |>
#   magrittr::set_names(c("pbmc3k1_1", "pbmc3k1_2", "pbmc3k1_3", "pbmc3k1_4")) 
#   
   dir("dev/CAQ_sce/", full.names = T) |> head(2)

  # Initialise pipeline characteristics
# file_list |> 
input_hpc |> 
  initialise_hpc(
    gene_nomenclature = "symbol",
    data_container_type = "sce_hdf5",
    store = "~/scratch/Census/temp5/",
    tier = c("tier_1","tier_1"),
    computing_resources = computing_resources,
    #debug_step ="empty_tbl_0cf8d597acd380df"
    # debug_step = "non_batch_variation_removal_S_1",


    # Default resourced 
  ) |> 
  
  # ONLY APPLICABLE TO SCE FOR NOW
  transform_assay(fx = file_list |> purrr::map(~identity), target_output = "sce_transformed") |> 
  
  # hpc_report(
  #   "empty_report", 
  #   rmd_path = paste0(system.file(package = "HPCell"), "/rmd/test.Rmd"), 
  #   empty_list = "empty_tbl" |> is_target(),
  #   sample_names = "sample_names" |> is_target()
  # ) |> 
  # 
  # 
  # hpc_iterate(
  #   target_output = "o",
  #   user_function = function(x, y){x |> dplyr::mutate(bla = y)},
  #   x = "data_object" |> is_target(),
  #   y = "works"
  # ) |>

  hpc_iterate(
    target_output = "foo",
    user_function = foo_function |> quote(),
    x = "data_object" |> is_target(),
    y = "works", 
    user_function_source_path = "dev/my_custom_script.R" |> here::here()
  ) |>
  
  # Remove empty outliers
  remove_empty_DropletUtils( target_input = "sce_transformed") |> 
  
  # Annotation
  annotate_cell_type(
    target_input = "data_object"
    # ,
    # azimuth_reference = readRDS("dev/reference_azimuth.rds")
  ) |> 
  
  # Remove dead cells
  remove_dead_scuttle(target_input = "data_object") |> 
  
  # Score cell cycle
  score_cell_cycle_seurat(target_input = "data_object") |> 
  
  # Remove doublets
  remove_doublets_scDblFinder(target_input = "data_object") |> 
  
  normalise_abundance_seurat_SCT(factors_to_regress = c(
    "subsets_Mito_percent",
    "subsets_Ribo_percent",
    "G2M.Score"
  ), 
  target_input = "data_object") |> 
  
  calculate_pseudobulk(group_by = "monaco_first.labels.fine", target_input = "data_object") |> 
  
  # test_differential_abundance(~ age_days + (1|collection_id), .abundance="counts") |> 
  test_differential_abundance(~ age_days, .abundance="counts", group_by_column = "monaco_first.labels.fine") |> 

  # For the moment only available for single cell
  get_single_cell(target_input = "data_object") 


## Report testing 
  
library(HPCell)
library(targets)
library(Seurat)
library(SeuratData)
library(crew)
library(crew.cluster)

bp<- scRNAseq::fetchDataset("baron-pancreas-2016", "2023-12-14", path="human") 
  
file.path <- ("~/HPCell/bp.rds")
split_values <- unique(bp$label)

# Create a list of SCE objects split by the column
sce_list <- lapply(split_values, function(value) {
  bp[, bp$label == value]
})


##### Testing args 
# empty_tbl <- tar_read(empty_tbl)
# data_object <- tar_read(data_object)
# alive_tbl <- tar_read(alive_tbl)
# sample_name <- tar_read(sample_names)
# cell_cycle_tbl <- tar_read(cell_cycle_tbl)
# annotation_tbl <- tar_read(annotation_tbl)
# doublet_tbl <- tar_read(doublet_tbl)
# data_object <- tar_read(data_object)
#########################################
library(HPCell)
library(targets)
library(Seurat)
library(SeuratData)
library(crew)
library(crew.cluster)

InstallData("ifnb")
ifnb <- UpdateSeuratObject(ifnb)
ifnb.list <- SplitObject(ifnb, split.by = "stim")
file_paths <- c("~/CTRL_seurat_tibble.rds", "~/STIM_seurat_tibble.rds")

#Subset 300 cells 
ctrl_subset <- subset(ifnb.list$CTRL, cells = sample(Cells(ifnb.list$CTRL), size = 300))
stim_subset <- subset(ifnb.list$STIM, cells = sample(Cells(ifnb.list$STIM), size = 300))

# Save the Seurat objects to the specified file paths
saveRDS(ctrl_subset, file_paths[1])  # Save CTRL object
saveRDS(stim_subset, file_paths[2])  # Save STIM object

input_hpc =
  file_paths |> 
  magrittr::set_names(c("CTRL", "STIM"))

input_hpc |> 
  initialise_hpc(
  gene_nomenclature = "symbol",
  data_container_type = "seurat_rds", 
  computing_resources = crew_controller_local(workers = 8),
  ) |> 
  remove_empty_DropletUtils() |>          # Remove empty outliers
  remove_dead_scuttle() |>                # Remove dead cells
  score_cell_cycle_seurat() |>            # Score cell cycle
  remove_doublets_scDblFinder() |>        # Remove doublets
  annotate_cell_type() |>                 # Annotation across SingleR and Seurat Azimuth
  normalise_abundance_seurat_SCT(factors_to_regress = c(
    "subsets_Mito_percent", 
    "subsets_Ribo_percent", 
    "G2M.Score"
  )) |> 

  hpc_report(
    "empty_report",
    rmd_path = system.file("rmd", "Empty_droptlet_report.qmd", package = "HPCell"),
    empty_tbl = "empty_tbl" |> is_target(),
    data_object = "data_object" |> is_target(),
    alive_tbl = "alive_tbl" |> is_target(),
    sample_name = "sample_names" |> is_target()
  ) 
# |>
  hpc_report(
    "doublet_report",
    rmd_path = system.file("rmd", "Doublet_identification_report.qmd", package = "HPCell"),
    data_object = "data_object" |> is_target(),
    doublet_tbl = "doublet_tbl" |> is_target(),
    annotation_tbl = "annotation_tbl" |> is_target(),
    sample_names = "sample_names" |> is_target()
  ) |>
  hpc_report(
    "Technical_variation_report",
    rmd_path = system.file("rmd", "technical_variation_report.qmd", package = "HPCell"),
    data_object = "data_object" |> is_target(),
    empty_tbl = "empty_tbl" |> is_target(),
    sample_name = "sample_names" |> is_target()
  ) |>
  hpc_report(
  "pseudo_bulk_report",
  rmd_path = system.file("rmd", "pseudobulk_analysis_report.qmd", package = "HPCell"),
  data_object = "data_object" |>  is_target(),
  empty_tbl = "empty_tbl" |> is_target(),
  alive_tbl = "alive_tbl" |> is_target(),
  cell_cycle_tbl = "cell_cycle_tbl" |> is_target(),
  annotation_tbl = "annotation_tbl" |> is_target(),
  doublet_tbl = "doublet_tbl" |> is_target(),
  sample_name = "sample_names" |> is_target()
  )





library(HPCell)
library(targets)
library(Seurat)
library(SeuratData)
library(crew)
library(crew.cluster)

file_paths <- file.path(getwd(), c("CTRL_seurat_tibble.rds", "STIM_seurat_tibble.rds"))
names(file_paths) <- c("CTRL", "STIM")

# Install and prepare IFNB demo dataset
SeuratData::InstallData("ifnb")
ifnb <- ifnb |> UpdateSeuratObject()
ifnb.list <- SplitObject(ifnb, split.by = "stim")

# Sample 300 cells per condition 
set.seed(42)
ctrl_subset <- subset(ifnb.list$CTRL, cells = sample(Cells(ifnb.list$CTRL), size = 300))
stim_subset <- subset(ifnb.list$STIM, cells = sample(Cells(ifnb.list$STIM), size = 300))

saveRDS(ctrl_subset, file_paths["CTRL"])
saveRDS(stim_subset, file_paths["STIM"])

input_hpc <- file_paths

input_hpc =
  file_paths |> 
  magrittr::set_names(c("CTRL", "STIM"))

input_hpc |> 
  initialise_hpc(
    gene_nomenclature = "symbol",
    data_container_type = "seurat_rds", 
    computing_resources = crew_controller_local(workers = 8),
  ) |> 
  remove_empty_DropletUtils() |>          # Remove empty outliers
  remove_dead_scuttle() |>                # Remove dead cells
  score_cell_cycle_seurat() |>            # Score cell cycle
  remove_doublets_scDblFinder() |>        # Remove doublets
  annotate_cell_type() |>                 # Annotation across SingleR and Seurat Azimuth
  normalise_abundance_seurat_SCT(factors_to_regress = c(
    "subsets_Mito_percent", 
    "subsets_Ribo_percent", 
    "G2M.Score"
  )) |> 
  hpc_report(
    "empty_report",
    rmd_path = system.file("rmd", "Empty_droptlet_report.qmd", package = "HPCell"),
    empty_tbl = "empty_tbl" |> is_target(),
    data_object = "data_object" |> is_target(),
    alive_tbl = "alive_tbl" |> is_target(),
    sample_name = "sample_names" |> is_target()
  ) |>
  hpc_report(
    "doublet_report",
    rmd_path = system.file("rmd", "Doublet_identification_report.qmd", package = "HPCell"),
    data_object = "data_object" |> is_target(),
    doublet_tbl = "doublet_tbl" |> is_target(),
    annotation_tbl = "annotation_tbl" |> is_target(),
    sample_names = "sample_names" |> is_target()
  ) |>
  hpc_report(
    "Technical_variation_report",
    rmd_path = system.file("rmd", "technical_variation_report.qmd", package = "HPCell"),
    data_object = "data_object" |> is_target(),
    empty_tbl = "empty_tbl" |> is_target(),
    sample_name = "sample_names" |> is_target()
  ) |>
  hpc_report(
    "pseudo_bulk_report",
    rmd_path = system.file("rmd", "pseudobulk_analysis_report.qmd", package = "HPCell"),
    data_object = "data_object" |>  is_target(),
    empty_tbl = "empty_tbl" |> is_target(),
    alive_tbl = "alive_tbl" |> is_target(),
    cell_cycle_tbl = "cell_cycle_tbl" |> is_target(),
    annotation_tbl = "annotation_tbl" |> is_target(),
    doublet_tbl = "doublet_tbl" |> is_target(),
    sample_name = "sample_names" |> is_target()
  )

