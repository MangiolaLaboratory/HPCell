library(zellkonverter)
library(dplyr)
library(stringr)
library(Seurat)
library(SummarizedExperiment)
library(Matrix)
library(purrr)
library(glue)

x = sample_names |> head(2)
f = functions |> head(2)
tr = feature_thresh |> head(2)
job::job({
  
  library(HPCell)
  
  x |>
    initialise_hpc(
      store = "~/scratch/test_sct_target_store",
      gene_nomenclature = "ensembl",
      data_container_type = "anndata",
      # tier = tiers, # WE DON"T NEED AS WE HAVE ELASTIC RESOURCES NOW
      computing_resources = list(
        
        crew.cluster::crew_controller_slurm(
          name = "elastic",
          workers = 300,
          tasks_max = 20,
          seconds_idle = 30,
          crashes_error = 10,
          options_cluster = crew.cluster::crew_options_slurm(
            #memory_gigabytes_required = c(20, 35, 50, 75, 100, 150), 
            #memory_gigabytes_required = c(90, 100, 120, 150, 200,240), 
            #memory_gigabytes_required = c(60, 80, 100, 150, 200), 
            memory_gigabytes_required = c(45, 60, 75, 100, 120, 150), 
            cpus_per_task = c(2, 2, 5, 10, 20), 
            time_minutes = c(60*24, 60*24, 60*24, 60*24, 60*24,60*24),
            verbose = T
          )
        )
        
      ),
      verbosity = "summary",
      update = "never", 
      #update = "thorough", 
      error = "continue",
      garbage_collection = 100, 
      workspace_on_error = TRUE
      
    ) |> 
    transform_assay(fx = f, target_output = "sce_transformed") |>
    
    # # Remove empty outliers based on RNA count threshold per cell
    remove_empty_threshold(target_input = "sce_transformed", RNA_feature_threshold = tr) |>
    
    # Annotation
    annotate_cell_type(target_input = "sce_transformed", azimuth_reference = "pbmcref") |> 
    
    # Cell type harmonisation
    celltype_consensus_constructor(target_input = "sce_transformed",
                                   target_output = "cell_type_concensus_tbl") |>
    
    # Alive identification
    remove_dead_scuttle(target_input = "sce_transformed", target_annotation = "cell_type_concensus_tbl",
                        group_by = "cell_type_unified_ensemble") |>
    
    # Doublets identification
    remove_doublets_scDblFinder(target_input = "sce_transformed") |>
    
    # SCT 
    normalise_abundance_seurat_SCT(target_input = "sce_transformed", factors_to_regress = c(
      "subsets_Mito_percent",
      "subsets_Ribo_percent")) |> 
    
    print()
})
tar_meta(store = my_store, starts_with("sct")) |>  dplyr::count(is.na(error))
tar_meta(store = my_store, starts_with("sct")) |> dplyr::filter(!is.na(error)) |> dplyr::count(error)
tar_meta(store = my_store, starts_with("sct")) |> dplyr::filter(!is.na(error)) |> saveRDS("~/scratch/cache_temp/hpcell_sct_mtx_error_target_tibble.rds")
tar_meta(store = my_store, starts_with("sct")) |> dplyr::filter(!is.na(error))  |>
  filter(error |> str_like("%incorrect number of dimensions%")) |> pull(name)

tar_meta(store = my_store, starts_with("sct")) |> dplyr::filter(!is.na(error))  |>
  filter(error |> str_like("%need at least 2 data points%")) |> pull(name)

tar_workspace(sct_matrix_98aebc9baff6e2a8, store = my_store)
sce_transformed <- sce_transformed |> left_join(empty_tbl, by = ".cell") |>
  dplyr::filter(!empty_droplet) 
sce_transformed =
  sce_transformed |>
  left_join(
    alive_tbl ,
    by=".cell"
  ) |> dplyr::filter(alive) 
sce_transformed =
  sce_transformed |>
  left_join(
    doublet_tbl ,
    by=".cell"
  ) |> dplyr::filter(scDblFinder.class != "doublet") 

sce_transformed |> dim()
debugonce(non_batch_variation_removal)
non_batch_variation_removal(sce_transformed, empty_tbl, alive_tbl, doublet_tbl, 
                            NULL, factors_to_regress = c("subsets_Mito_percent",
                                                         "subsets_Ribo_percent"),
                            external_path = "~/scratch/cache_temp", container_type = data_container_type)

# sample id b4683e7aff5e496be9e1782a1fc21aa8 returns 2 variable genes in SCTransform based on NB model, which made SCTransform failed.
# Inspecting the reason why this is an issue

sce_transformed |> assay("X") |> as.numeric() |> summary()
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0         0         0         2         0 485165194 

# Check counts from raw sample
sliced_sample_tbl <- readRDS("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/sliced_sample_tbl_2024_Jul.rds")
sample_anndata <- sliced_sample_tbl |> filter(sample_2 == "b4683e7aff5e496be9e1782a1fc21aa8") |> pull(file_name) |> 
  readH5AD(.x, reader = "R", use_hdf5 = T)
sample_anndata  |> assay("X") |> as.numeric() |> summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    0.00    0.08    0.00 3684.51

raw_counts = sample_anndata  |> assay("X")
raw_counts <- raw_counts|>  as("dgCMatrix")


var <- DelayedArray::rowVars(raw_counts)
var |> summary()

sample_anndata_filtered <- sample_anndata |> left_join(empty_tbl, by = ".cell") |>
  dplyr::filter(!empty_droplet)  |>
  left_join(
    alive_tbl ,
    by=".cell"
  ) |> dplyr::filter(alive) |>
  left_join(
    doublet_tbl ,
    by=".cell"
  ) |> dplyr::filter(scDblFinder.class != "doublet") 
raw_counts_filtered <- sample_anndata_filtered|> assay("X")|>  as("dgCMatrix")
genes_amean <- rowMeans(raw_counts_filtered)
genes_var <- row_var(raw_counts_filtered)
overdispersion_factor <- genes_var - genes_amean
# pseudo give all genes
genes_step1 <- rownames(raw_counts_filtered)
overdispersion_factor_step1 <- overdispersion_factor[genes_step1]
is_overdispersed <- (overdispersion_factor_step1 > 0)
is_overdispersed |> sum()


sce_transformed_filtered <- sce_transformed |> left_join(empty_tbl, by = ".cell") |>
  dplyr::filter(!empty_droplet)  |>
  left_join(
    alive_tbl ,
    by=".cell"
  ) |> dplyr::filter(alive) |>
  left_join(
    doublet_tbl ,
    by=".cell"
  ) |> dplyr::filter(scDblFinder.class != "doublet") 
var2 <- assay(sce_transformed_filtered, "X")|>  as("dgCMatrix")

genes_amean <- rowMeans(var2)
genes_var <- row_var(var2)
overdispersion_factor <- genes_var - genes_amean
# pseudo give all genes
genes_step1 <- rownames(var2)
overdispersion_factor_step1 <- overdispersion_factor[genes_step1]
is_overdispersed <- (overdispersion_factor_step1 > 0)
is_overdispersed |> sum()

# SCE transformation might destroy good things
# Inspect the distribution of samples those having overdispersion gene issues
sample_target_tbl = tar_meta(starts_with("sct_"), store = "/vast/scratch/users/shen.m/cellNexus_target_store_2024_Jul") |> 
  filter(!is.na(error))  |> filter(error |> str_detect("subscript out of bounds|incorrect number of dimensions|need at least 2 data points")) |>
  pull(name) |> 
  purrr::map( ~ {
    
    e <- new.env(parent = globalenv())
    
    do.call(
      targets::tar_workspace,
      list(
        name = as.name(.x),     # <-- convert string -> symbol
        envir = e,
        packages = TRUE,
        source = TRUE,
        store = "/vast/scratch/users/shen.m/cellNexus_target_store_2024_Jul"
      )
    )
    
    sce <- e$sce_transformed
    sampleId <- sce |> distinct(sample_id) |> pull()
    
    tibble(target_name = .x, sample_id = sampleId)
  }, .progress = T)

sample_target_tbl = sample_target_tbl |> bind_rows() |> mutate(sample_id = sample_id |> as.character())
sample_target_tbl = sample_target_tbl |> mutate(path = file.path("/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/2024-07-01/", paste0(sample_id, ".h5ad")))
#sample_target_tbl |> write_parquet("~/scratch/samples.parquet")

# Get summary table in parallel
summary_store = "/vast/scratch/users/shen.m/sct_failed_samples_raw_counts_summary_target_store"
tar_script({
  library(dplyr)
  library(SummarizedExperiment)
  library(zellkonverter)
  library(crew)
  library(crew.cluster)
  library(duckdb)
  
  
  computing_resources = crew_controller_slurm(
    name = "elastic",
    workers = 300,
    tasks_max = 20,
    seconds_idle = 30,
    crashes_error = 10,
    options_cluster = crew_options_slurm(
      #memory_gigabytes_required = c(30, 45, 45, 90, 200, 1000, 1500), 
      memory_gigabytes_required = c(90, 130, 170, 200, 1000, 1500), 
      cpus_per_task = c(2, 2, 2, 2, 2, 2), 
      time_minutes = c(30*4, 30*4, 30*4, 60*4, 60*24, 60*24),
      verbose = T
    ))
  
  tar_option_set(
    memory = "transient",
    garbage_collection = TRUE,
    storage = "worker",
    retrieval = "worker",
    format = "qs",
    cue = tar_cue(mode = "never"),
    error = "continue",
    controller = computing_resources
  )
  
  get_sample_summary_stats <- function(files) {
    sce = readH5AD(files, reader = "R", use_hdf5 = T)
    
    if (ncol(sce) == 0) return(NULL)
    
    assay_name = sce@assays |> names() |> magrittr::extract(1)
    counts_mat = sce |> assay(assay_name) |> as.matrix()
    counts_vec = as.numeric(counts_mat)
    sample_id = basename(files)
    
    # Perform checks
    min_val = min(counts_vec, na.rm = TRUE)
    max_val = max(counts_vec, na.rm = TRUE)
    has_negative = min_val < 0
    max_gt_10 = max_val > 10
    
    # Check if all integers
    tol = 1e-5
    all_integer = all(counts_vec == floor(counts_vec), na.rm = TRUE)
    
    #should've excluded integer, but I will fix the output dataframe
    #has_floating = all(abs(counts_vec - round(counts_vec)) < tol, na.rm = TRUE)
    has_floating = !all_integer && all(abs(counts_vec - round(counts_vec)) < tol, na.rm = TRUE)
    
    tbl = tibble::tibble(
      sample_id = sample_id,
      min_val = min_val,
      max_val = max_val,
      has_negative = has_negative,
      max_gt_10 = max_gt_10,
      all_integer = all_integer,
      has_floating = has_floating,
      n_cells = ncol(counts_mat),
      n_genes = nrow(counts_mat)
    )
    
    
    tbl
  }
  
  list(
    tar_target(
      files,
      arrow::read_parquet("~/scratch/samples.parquet") |> pull(path),
      deployment = "main"
    ),
    
    # Get raw counts matrix with sample_id
    tar_target(
      sample_summary_df,
      get_sample_summary_stats(files),
      pattern = map(files),
      iteration = "list"
    )
  )
}, ask = FALSE,  script = glue("{summary_store}/_targets.R"))


job::job({
  
  tar_make(
    # callr_function = NULL,
    reporter = "summary",
    script = glue("{summary_store}/_targets.R"),
    store = glue("{summary_store}/_targets")
  )
  
})

sample_summary_df =  tar_read(sample_summary_df, store = glue("{summary_store}/_targets")) |> bind_rows() |>
  mutate(max_gt_20 = ifelse(max_val > 20, TRUE, FALSE)) |> bind_rows()

# Distribution PDF
plot_raw_counts_hist_pdf <- function(
    samples_to_plot,
    h5ad_dir,
    pdf_file,
    width = 8,
    height = 8,
    nrow = 3,
    ncol = 3,
    max_cells = 5e3,
    transform = c("log1p", "none"),
    hist_ylim = c(0, 1e5)
) {
  
  transform <- match.arg(transform)
  
  pdf(pdf_file, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  
  par(mfrow = c(nrow, ncol), mar = c(3, 3, 2, 1))
  
  samples_to_plot |>
    dplyr::pull(sample_id) |>
    purrr::map(
      ~ {
        sce <- zellkonverter::readH5AD(
          file = file.path(h5ad_dir, .x),
          reader = "R",
          use_hdf5 = T
        )
        
        if (ncol(sce) == 0) return(NULL)
        
        if (ncol(sce) > max_cells) {
          sce <- dplyr::sample_n(sce, max_cells)
        }
        
        assay_name <- names(SummarizedExperiment::assays(sce))[1]
        
        dataset_id <- sce |>
          dplyr::distinct(dataset_id) |>
          dplyr::pull()
        
        x <- SummarizedExperiment::assay(sce, assay_name) |>
          as.numeric()
        
        if (transform == "log1p") {
          x <- log1p(x)
        }
        
        hist(
          x,
          main = dataset_id,
          ylim = hist_ylim
        )
        
        invisible(NULL)
      },
      .progress = TRUE
    )
  
  par(mfrow = c(1, 1))
}

# CASE 1
samples_to_plot <- sample_summary_df |> filter(!has_negative, !max_gt_20, !all_integer, !has_floating) |> 
  arrange(n_cells) |>
  slice_sample(n = 5) |>
  ungroup()
plot_raw_counts_hist_pdf(samples_to_plot, 
                         h5ad_dir = "/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/2024-07-01/",
                         pdf_file = "~/samples_case1.pdf", 
                         transform = "none",
                         hist_ylim = c(0, 1e4))


# CASE 2
samples_to_plot <- sample_summary_df |> filter(!has_negative, max_gt_20, !all_integer, !has_floating) |> 
  arrange(n_cells) |>
  slice_sample(n = 5) |>
  ungroup()
plot_raw_counts_hist_pdf(samples_to_plot, 
                         h5ad_dir = "/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/2024-07-01/",
                         pdf_file = "~/samples_case2.pdf", 
                         transform = "log1p",
                         hist_ylim = c(0, 1e4))



# CASE 3
samples_to_plot <- sample_summary_df |> filter(!has_negative, max_gt_20, all_integer, !has_floating) |> 
  arrange(n_cells) |>
  slice_sample(n = 5) |>
  ungroup()
plot_raw_counts_hist_pdf(samples_to_plot, 
                         h5ad_dir = "/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/2024-07-01/",
                         pdf_file = "~/samples_case3.pdf", 
                         transform = "log1p",
                         hist_ylim = c(0, 1e4))

# CASE 4
samples_to_plot <- sample_summary_df |> filter(has_negative, !max_gt_20, !all_integer, !has_floating) |> 
  arrange(n_cells) |>
  slice_sample(n = 5) |>
  ungroup()
plot_raw_counts_hist_pdf(samples_to_plot, 
                         h5ad_dir = "/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/2024-07-01/",
                         pdf_file = "~/samples_case4.pdf", 
                         transform = "none",
                         hist_ylim = c(0, 1e4))

# CASE 5
samples_to_plot <- sample_summary_df |> filter(has_negative, max_gt_20, !all_integer, !has_floating) |> 
  arrange(n_cells) |>
  slice_sample(n = 5) |>
  ungroup()
plot_raw_counts_hist_pdf(samples_to_plot, 
                         h5ad_dir = "/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/2024-07-01/",
                         pdf_file = "~/samples_case5.pdf", 
                         transform = "none",
                         hist_ylim = c(0, 1e4))

# Decision:
# Samples failed sctransform due to wrong transformation label are re-evaluated in ~/git_control/HPCell/dev/cellnexus-2024-scripts/samples_reannotate_transformation_for_fail_sct_cellnexus2024.R,

