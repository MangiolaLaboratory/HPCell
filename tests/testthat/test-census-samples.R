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
store = "~/scratch/Census/census_reanalysis/census-run-samples"
files <- dir(glue("{directory}"), full.names = T)
# results <- purrr::map_dfr(files, function(file_path) {
#   data <- zellkonverter::readH5AD(file_path, use_hdf5 = TRUE, reader = "R", verbose = TRUE)
# 
#   cell_number <- length(colnames(data))
# 
# 
#   file_size <- file.info(file_path)$size / 1024^3
# 
#   # nFeature_threshold
# 
#   tibble(file_name = file_path,
#          cell_number = cell_number,
#          file_size = file_size)
# })

#results <- readRDS(glue("{store}/sample_tiers.rds"))
results <- tar_read(results, store = glue("{store}/sample_tiers_dataframe_targets"))

tiers_dataframe <- results |>
  mutate(file_name = file_name,
         file_size = round(file_size, 3),
         tier = case_when(cell_number < 6000 ~ "tier_1",
                          cell_number > 6000 & cell_number <= 10000 ~ "tier_2",
                          cell_number > 10000 & cell_number < 20000 ~ "tier_3",
                          cell_number > 20000 & cell_number < 40000 ~ "tier_4",
                          cell_number > 40000 ~ "tier_5"),
         set_names = basename(file_name) |> stringr::str_remove("\\.h5ad$"))

result_directory = "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024"
samples <- read_parquet("~/cellxgene_curated/census_samples/census_samples_to_download_groups.parquet")
sample_meta <- tar_read(metadata_dataset_id_common_sample_columns, store = glue("{result_directory}/_targets"))
samples = samples |> left_join(get_metadata() |> select(dataset_id, contains("norm")) |>
                                 distinct() |> filter(!is.na(x_normalization)) |>
                                 as_tibble(), by = "dataset_id")


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
    ~ ( function(data) {
      assay_name <- data@assays |> names() |> magrittr::extract2(1)
      counts <- assay(data, assay_name)
      density_est <- density(counts |> HPCell:::get_count_per_gene_df() |> pull(counts) )
      mode_value <- density_est$x[which.max(density_est$y)]
      if (mode_value < 0 )  counts <- counts + abs(mode_value)
      
      # Scale max counts to 20 to avoid any downstream failure
      if ((.x == "exp") && (max(counts) > 20)) {
        scale_factor = 20 / max(counts)
        counts <- counts * scale_factor}
      
      counts <- transform_method(counts)
      # round counts to avoid potential substraction error due to different digits print out 
      counts <- counts |> round(5)
      majority_gene_counts = names(which.max(table(as.vector(counts)))) |> as.numeric()
      if (majority_gene_counts != 0) {
        counts <- counts - majority_gene_counts
      }
      
      # Avoid downstream failures negative counts
      if((counts[,1:min(10000, ncol(counts))] |> min()) < 0)
        counts[counts < 0] <- 0
      
      col_sums <- colSums(counts)
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


files <- results |> mutate(sample_2 = basename(file_name) |> stringr::str_remove("\\.h5ad$")) |> 
  left_join(df, by = "sample_2") |> left_join(tiers_dataframe, by = c("sample_2"= "set_names")) |> 
  select(-file_name.y, -cell_number.y, -file_size.y) |> rename(file_name = file_name.x,
                                                               cell_number = cell_number.x,
                                                               file_size = file_size.x)


file_list = files |> head(100)

file_list |> pull(file_name) |>
  initialise_hpc(
    gene_nomenclature = "ensembl",
    data_container_type = "anndata",
    store = "~/scratch/Census/census_reanalysis/census-run-samples/50samples_pilot/",
    #debug_step = "empty_tbl_69e0f0d88e44d787",
    tier = file_list |> pull(tier),
    computing_resources = list(

      crew_controller_slurm(
        name = "tier_1",
        slurm_memory_gigabytes_per_cpu = 20,
        slurm_cpus_per_task = 1,
        workers = 50,
        tasks_max = 5,
        verbose = T
      ),
      crew_controller_slurm(
        name = "tier_2",
        slurm_memory_gigabytes_per_cpu = 30,
        slurm_cpus_per_task = 1,
        workers = 50,
        tasks_max = 5,
        verbose = T
      ),
      crew_controller_slurm(
        name = "tier_3",
        slurm_memory_gigabytes_per_cpu = 50,
        slurm_cpus_per_task = 1,
        workers = 50,
        tasks_max = 5,
        verbose = T
      ),
      crew_controller_slurm(
        name = "tier_4",
        slurm_memory_gigabytes_per_cpu = 100,
        slurm_cpus_per_task = 1,
        workers = 50,
        tasks_max = 5,
        verbose = T
      ),
      crew_controller_slurm(
        name = "tier_5",
        slurm_memory_gigabytes_per_cpu = 200,
        slurm_cpus_per_task = 1,
        workers = 50,
        tasks_max = 5,
        verbose = T
      )
      )
    
  ) |> 
  #tranform_assay(fx =  purrr::map(1:20, ~identity), target_output = "sce_transformed") |> 
  tranform_assay(fx = file_list |>
                   pull(transformation_function),
                 target_output = "sce_transformed") |>
  
  # Remove empty outliers based on RNA count threshold per cell
  remove_empty_threshold(target_input = "sce_transformed", RNA_feature_threshold = 200 ) |>
  
  # Remove empty outliers
  #remove_empty_DropletUtils(target_input = "sce_transformed") |>
  
  # Remove dead cells
  remove_dead_scuttle(target_input = "sce_transformed") |>
  
  # Score cell cycle
  score_cell_cycle_seurat(target_input = "sce_transformed") |>
  
  # Remove doublets
  remove_doublets_scDblFinder(target_input = "sce_transformed") |>
  
  # Annotation
  annotate_cell_type(target_input = "sce_transformed", azimuth_reference = "pbmcref") |>
  
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

  

