# test census-20-samples --------------------------------------------------
library(glue)
library(dplyr)
library(arrow)
library(SingleCellExperiment)
library(tibble)
library(SummarizedExperiment)
library(CuratedAtlasQueryR)
library(glue)
library(purrr)
library(zellkonverter)
library(tidyr)
library(ggplot2)
library(plotly)
library(targets)
library(stringr)
directory = "~/cellxgene_curated/census_samples/anndata"
store = "~/scratch/Census/census_reanalysis/census-run-20-samples"
files <- dir(glue("{directory}"), full.names = T) |> head(20)
# results <- purrr::map_dfr(files, function(file_path) {
#   data <- zellkonverter::readH5AD(file_path, use_hdf5 = TRUE, reader = "R", verbose = TRUE)
#   
#   cell_number <- length(colnames(data))
#   
#   
#   file_size <- file.info(file_path)$size / 1073741824
#   
#   tibble(file_name = file_path,
#          cell_number = cell_number,
#          file_size = file_size)
# })

#results |> saveRDS(glue("{store}/sample_tiers.rds"))
results <- readRDS(glue("{store}/sample_tiers.rds"))

tiers_dataframe <- results |>
  mutate(file_name = file_name,
         file_size = round(file_size, 3),
         tier = ifelse(cell_number > 6000, "tier_2", "tier_1"),
         set_names = basename(file_name) |> stringr::str_remove("\\.h5ad$"))

result_directory = "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024"
samples <- read_parquet("~/cellxgene_curated/census_samples/census_samples_to_download_groups.parquet")
sample_meta <- tar_read(metadata_dataset_id_common_sample_columns, store = glue("{result_directory}/_targets"))
samples = samples |> left_join(get_metadata() |> select(dataset_id, contains("norm")) |>
                                 distinct() |> filter(!is.na(x_normalization)) |>
                                 as_tibble(), by = "dataset_id")
# df <- samples |> left_join(sample_meta, by = "dataset_id") |> distinct(dataset_id, sample_2, x_normalization, x_approximate_distribution) |>
#   mutate(current_normalisation_method = case_when(str_like(x_normalization, "C%") ~ "log",
#                                                   x_normalization == "none" ~ "log",
#                                                   x_normalization == "normalized" ~ "log",
#                                                   is.na(x_normalization) & is.na(x_approximate_distribution) ~ "log",
#                                                   is.na(x_normalization) & x_approximate_distribution == "NORMAL" ~ "NORMAL",
#                                                   is.na(x_normalization) & x_approximate_distribution == "COUNT" ~ "COUNT",
#                                                   str_like(x_normalization, "%canpy%") ~ "log1p",
#                                                   
#                                                   TRUE ~ x_normalization)) |>
#   mutate(method_to_apply = case_when(current_normalisation_method %in% c("log","LogNormalization","LogNormalize","log-normalization") ~ "exp(counts)",
#                                      is.na(x_normalization) & is.na(x_approximate_distribution) ~ "exp(counts)",
#                                      str_like(current_normalisation_method, "Counts%") ~ "exp(counts)",
#                                      str_like(current_normalisation_method, "%log2%") ~ "2^counts",
#                                      current_normalisation_method %in% c("log1p", "log1p, base e", "Scanpy", 
#                                                                          "scanpy.api.pp.normalize_per_cell method, scaling factor 10000") ~ "expm1(counts)",
#                                      current_normalisation_method == "log1p, base 2" ~ "2^counts - 1",
#                                      current_normalisation_method == "NORMAL" ~ "exp(counts)",
#                                      current_normalisation_method == "COUNT" ~ "no method applied"
#   )
#   ) |>
#   mutate(comment = case_when(str_like(x_normalization, "Counts%")  ~ "a checkpoint for max value of Assay must <= 50",
#                              is.na(x_normalization) & is.na(x_approximate_distribution) ~ "round negative value to 0",
#                              x_normalization == "normalized" ~ "round negative value to 0"
#                              
#   ))

df <- samples |> left_join(sample_meta, by = "dataset_id") |> distinct(dataset_id, sample_2, x_normalization, x_approximate_distribution) |>
  mutate(transform_method = case_when(str_like(x_normalization, "C%") ~ "log",
                                      x_normalization == "none" ~ "log",
                                      x_normalization == "normalized" ~ "log",
                                      is.na(x_normalization) & is.na(x_approximate_distribution) ~ "log",
                                      is.na(x_normalization) & x_approximate_distribution == "NORMAL" ~ "NORMAL",
                                      is.na(x_normalization) & x_approximate_distribution == "COUNT" ~ "COUNT",
                                      str_like(x_normalization, "%canpy%") ~ "log1p",
                                      TRUE ~ x_normalization)) |>
  
  mutate(method_to_apply =  case_when(transform_method %in% c("log","LogNormalization","LogNormalize","log-normalization") ~ "exp",
                                      is.na(x_normalization) & is.na(x_approximate_distribution) ~ "exp",
                                      str_like(transform_method, "Counts%") ~ "exp",
                                      str_like(transform_method, "%log2%") ~ "exp",
                                      transform_method %in% c("log1p", "log1p, base e", "Scanpy", 
                                                              "scanpy.api.pp.normalize_per_cell method, scaling factor 10000") ~ "expm1",
                                      transform_method == "log1p, base 2" ~ "expm1",
                                      transform_method == "NORMAL" ~ "exp",
                                      transform_method == "COUNT" ~ "identity"
  ) ) |>
  mutate(comment = case_when(str_like(x_normalization, "Counts%")  ~ "a checkpoint for max value of Assay must <= 50",
                             is.na(x_normalization) & is.na(x_approximate_distribution) ~ "round negative value to 0",
                             x_normalization == "normalized" ~ "round negative value to 0"
  )) |> 
  mutate(transformation_function = map(
    method_to_apply,
    ~ (function(data) {
      assay_name <- data@assays |> names() |> magrittr::extract2(1)
      counts <- assay(data, assay_name)
      density_est <- density(counts |> HPCell:::get_count_per_gene_df() |> pull(counts) )
      mode_value <- density_est$x[which.max(density_est$y)]
      if (mode_value < 0 )  counts <- counts + abs(mode_value)
      
      # Scale max counts to 20 to avoid any downstream failure
      if ((.x == "exp") && (max(counts) > 20)){
        scale_factor = 20/max(counts)
        counts = counts * scale_factor
        # Apply the transformation
        counts <- transform_method(counts)
        
        # Avoid majority of genes after transformation are 1, so substract by 1
        majority_gene_counts = names(which.max(table(as.vector(counts)))) |> as.numeric()
        
        #substract all counts by the majority_gene_counts if majority_gene_counts is not 0.
        if (majority_gene_counts != transform_method(scale_factor)) counts <- counts - majority_gene_counts
      }
      

      # Avoid downstream failures negative counts
      if((counts[,1:min(10000, ncol(counts))] |> min()) < 0)
        counts[counts < 0] <- 0
      
      col_sums <- colSums(counts)
      # Cap large values
      if(max(col_sums) > 1e100) {
        temp <- counts[, sample(1:ncol(data), size = 10000, replace = TRUE), drop = TRUE]
        q <- quantile(temp[temp > 0], 0.9)
        counts[counts > 1e100] <- q
      }
      # Drop all zero cells
      data <- data[, col_sums > 0]
      
      # Avoid downstream binding error
      rowData(data) = NULL
      
      # Assign counts back to data
      assay(data, assay_name) <- counts
      
      data
      
    }) |>
      # Meta programming, replacing the transformation programmatically
      substitute( env = list(transform_method = as.name(.x))) |>
      # Evaluate back to a working function
      eval()
  ))


files <- results |> mutate(sample_2 = basename(file_name) |> tools::file_path_sans_ext()) |> 
  left_join(df, by = "sample_2") |> left_join(tiers_dataframe, by = c("sample_2"= "set_names")) |> 
  select(-file_name.y, -cell_number.y, -file_size.y) |> rename(file_name = file_name.x,
                                                               cell_number = cell_number.x,
                                                               file_size = file_size.x)
files |> head(20) |> pull(file_name) |>
  initialise_hpc(
    gene_nomenclature = "ensembl",
    data_container_type = "anndata",
    store = "~/scratch/Census/census_reanalysis/census-run-20-samples/20samples/",
    #debug_step = "cell_cycle_tbl_tier_1_e0606bff6b731c54",  # problematic sample
    #debug_step = "cell_cycle_tbl_tier_1_ea6705b8ecef4374",
    tier =  files |> head(20)  |> pull(tier),
    #tier = "tier_1",
    #computing_resources = crew_controller_local(workers = 10) #resource_tuned_slurm
    computing_resources = list(

      crew_controller_slurm(
        name = "tier_1",
        slurm_memory_gigabytes_per_cpu = 15,
        slurm_cpus_per_task = 1,
        workers = 50,
        tasks_max = 5,
        verbose = T,
        slurm_time_minutes = 120
      ),
      crew_controller_slurm(
        name = "tier_2",
        slurm_memory_gigabytes_per_cpu = 30,
        slurm_cpus_per_task = 1,
        workers = 50,
        tasks_max = 5,
        verbose = T,
        slurm_time_minutes = 120
      ))
    
  ) |> 
  # this does not tested whether identity being successfully assigned or unassigned, because identity remains the original function
  #tranform_assay(fx =  purrr::map(1:20, ~identity), target_output = "sce_transformed") |> 
  # tranform_assay(fx =  files |> filter(file_name == "/home/users/allstaff/shen.m/cellxgene_curated/census_samples/anndata/0024f909bf540734c021854ee1c758ca.h5ad") |> 
  #                  pull(transformation_function) |> _[[1]], 
  #                  target_output = "sce_transformed") |>
  tranform_assay(fx =   files |> head(20) |>
                   pull(transformation_function),
                 target_output = "sce_transformed") |>
  
  # Remove empty outliers based on RNA count threshold per cell
  remove_empty_threshold(target_input = "sce_transformed", RNA_count_threshold = 100 ) |>
  
  # Remove empty outliers
  # remove_empty_DropletUtils(target_input = "sce_transformed") |>
  
  # Remove dead cells
  remove_dead_scuttle(target_input = "sce_transformed") |>
  
  # Score cell cycle
  score_cell_cycle_seurat(target_input = "sce_transformed") |>
  
  # Remove doublets
  remove_doublets_scDblFinder(target_input = "sce_transformed") |>
  
  # Annotation
  annotate_cell_type(target_input = "sce_transformed") |>
  
  normalise_abundance_seurat_SCT(factors_to_regress = c(
    "subsets_Mito_percent",
    "subsets_Ribo_percent",
    "G2M.Score"
  ), target_input = "sce_transformed")
  
   # calculate_pseudobulk(group_by = "monaco_first.labels.fine", target_input = "sce_transformed") |>
   # 
   # # test_differential_abundance(~ age_days + (1|collection_id), .abundance="counts") |>
   # # #test_differential_abundance(~ age_days, .abundance="counts")
   # #
   # # # For the moment only available for single cell
   # get_single_cell(target_input = "sce_transformed")


# # Test substitute eval function
#   .x = "exp"
#   fx <- (function(x){ bla(x) }) |> substitute( env = list(bla = as.name(.x))) |> eval()
#   fx(2)
#   
#   
#   get_count_per_gene_df <- function(counts) {
#     counts_tidy <- counts |> as.data.frame() |> tibble::rownames_to_column(var = "features") |>
#       as_tibble() |> pivot_longer(!features, names_to = "cells",
#                                   values_to = "counts")
#     counts_tidy
#   }
#   
#   .x = "exp"
#   my_function <- (function(data) {
#     assay_name <- data@assays |> names() |> magrittr::extract2(1)
#     counts <- assay(data, assay_name)
#     density_est <- density(counts |> get_count_per_gene_df() |> pull(counts) )
#     mode_value <- density_est$x[which.max(density_est$y)]
#     if (mode_value < 0 )  counts <- counts + abs(mode_value)
#     # Apply the transformation
#     counts <- transform_method(counts)
#     
#     # Avoid downstream failures
#     if((counts[,1:min(10000, ncol(counts))] |> min()) < 0)
#       counts[counts < 0] <- 0
#     
#     # Avoid majority of genes after transformation are 1, so substract by 1
#     counts = names(which.max(table(as.vector(counts))))
#     
#     if (counts == "1") counts <- counts - 1
#     col_sums <- colSums(counts)
#     # Cap large values
#     if(max(col_sums) > 1e100) {
#       temp <- counts[, sample(1:ncol(data), size = 10000, replace = TRUE), drop = TRUE]
#       q <- quantile(temp[temp > 0], 0.9)
#       counts[counts > 1e100] <- q
#     }
#     # Drop all zero cells
#     data <- data[, col_sums > 0]
#     
#     # Avoid downstream binding error
#     rowData(data) = NULL
#     
#     # Assign counts back to data
#     assay(data, assay_name) <- counts
#     
#     data
#     
#   }) |>
#     # Meta programming, replacing the transformation programmatically
#     substitute( env = list(transform_method = as.name(.x))) |>
#     # Evaluate back to a working function
#     eval()
#   
#   original_data = readH5AD("/home/users/allstaff/shen.m/cellxgene_curated/census_samples/anndata/0024f909bf540734c021854ee1c758ca.h5ad", reader = "R", use_hdf5 = TRUE)
#   original_data |> assay("X") |> as.numeric() |> summary()
#   transform_data <- my_function(original_data)
#   transform_data |> assay("X") |> as.numeric() |> summary()
#   
#   
#   transform_data |> assay("X") |> get_count_per_gene_df() |>  ggplot(aes(counts)) +geom_density()
#   
  
  

