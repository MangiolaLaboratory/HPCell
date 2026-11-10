# This script reannotate transformation for samples that failed SCTransform (v2) in 
#    HPCell, due to fewer overdispersion genes were diagnosed.
# Debug overdispersion genes script is ~/git_control/HPCell/dev/cellnexus-2024-scripts/run_sct_for_cellNexus_test_pipeline_debug_failing_sct.R

library(zellkonverter)
library(dplyr)
library(stringr)
library(Seurat)
library(SummarizedExperiment)
library(Matrix)
library(purrr)
library(targets)
library(arrow)

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

sliced_sample_tbl <- readRDS("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/sliced_sample_tbl_2024_Jul.rds")

sample_summary_df = sample_summary_df |> mutate(sample_id = str_remove(sample_id, ".h5ad"))
impute_x_approximate_distribution <- function(df) {
  df |> mutate(
    inferred_distribution = case_when(

      # 1) No negatives, no large values, no integers, no floating
      !has_negative & !max_gt_20 & !all_integer & !has_floating ~ "log1p",
      
      # 2) No negatives, has large values
      !has_negative &  max_gt_20 & !all_integer & !has_floating  ~ "raw",
      
      # 3) No negatives, large values, all integer, has floating
       !has_negative &  max_gt_20 & all_integer & !has_floating ~ "raw",
      
      # 4) Has negatives, no large values, no integer, no floating
      has_negative & !max_gt_20 & !all_integer & !has_floating ~ "raw",
      
      # 5) Has negatives and large values
      has_negative &  max_gt_20 & !all_integer & !has_floating ~ "raw"
    )
  )
} 

sample_summary_df = sample_summary_df |> impute_x_approximate_distribution() |> 
  # Inverse distribution
  mutate(method_to_apply = case_when(inferred_distribution == "log" ~ "exp",
                                     inferred_distribution == "log1p" ~ "expm1",
                                     inferred_distribution == "raw" ~ "identity"))
  
sliced_sample_tbl = sliced_sample_tbl |> left_join(sample_summary_df |> mutate(sample_id = str_remove(sample_id, ".h5ad")) |>
                                 select(sample_id, method_to_apply), by = c("sample_2" =  "sample_id")) |>
  mutate(method_to_apply = case_when(is.na(method_to_apply.y) ~ method_to_apply.x,
                                     !is.na(method_to_apply.y) ~ method_to_apply.y))

sliced_sample_tbl |> saveRDS("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/sliced_sample_tbl_2024_Jul.rds")

# TEST RUNNING THESE SAMPLES IN HPCEll
saved <- sample_summary_df |> write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/samples_updated_transformation_method_to_fix_sct_2024_Jul.parquet")
updated_transformation_sliced_samples <- sliced_sample_tbl |> filter(sample_2 %in% (
  sample_summary_df |> pull(sample_id)
))

x <-
  updated_transformation_sliced_samples |> 
  pull(file_name) |> 
  str_replace("/home/users/allstaff/shen.m/scratch", "/vast/scratch/users/shen.m") |> 
  set_names(updated_transformation_sliced_samples |> pull(sample_2))
f  = updated_transformation_sliced_samples |> pull(method_to_apply)
tr = updated_transformation_sliced_samples |> pull(feature_thresh)

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
            memory_gigabytes_required = c(45, 60, 75, 100, 120, 150), 
            cpus_per_task = c(2, 2, 5, 10, 20), 
            time_minutes = c(60*24, 60*24, 60*24, 60*24, 60*24,60*24),
            verbose = T
          )
        )
        
      ),
      verbosity = "summary",
      update = "thorough", 
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

tar_meta( starts_with("sct_"), store = "~/scratch/test_sct_target_store" ) |> 
  dplyr::count(!is.na(error))

tar_workspace(sct_matrix_03dd1be4995c69ac,store = "~/scratch/test_sct_target_store" )
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
genes_var <- DelayedArray::rowVars(var2)
overdispersion_factor <- genes_var - genes_amean

hist(overdispersion_factor)
var2 |> as.numeric() |> hist()

# pseudo give all genes
genes_step1 <- rownames(var2)
overdispersion_factor_step1 <- overdispersion_factor[genes_step1]
is_overdispersed <- (overdispersion_factor_step1 > 0)
is_overdispersed |> sum() # WHY STILL 0 overdispersion gene

sample_anndata <- sliced_sample_tbl |> filter(sample_2 == "e560e83e58be239de24d583f692a81a1") |> pull(file_name) |> 
  readH5AD(.x, reader = "R", use_hdf5 = T)
sample_anndata  |> assay("X") |> as.numeric() |> hist()


# 0 overdispersion genes persists for some samples. Check their abundance
samples_to_inspect = tar_meta( starts_with("sct_"), store = "~/scratch/test_sct_target_store" ) |> 
  dplyr::filter(!is.na(error)) |>
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
        store ="~/scratch/test_sct_target_store"
      )
    )
    
    sce <- e$sce_transformed
    sampleId <- sce |> distinct(sample_id) |> pull()
    
    tibble(target_name = .x, sample_id = sampleId)
  }, .progress = T)

samples_to_inspect = samples_to_inspect |> bind_rows()
samples_to_inspect


sliced_sample_tbl <- readRDS("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/sliced_sample_tbl_2024_Jul.rds")
sliced_sample_tbl |> filter(sample_2 %in% (samples_to_inspect |> pull(sample_id))) |>
  filter(cell_number >1) |> head(1) |> pull(file_name) |> 
  readH5AD(.x, reader = "R", use_hdf5 = T)|> assay("X") |>
  as.numeric() |> hist()


# iteration
files_to_plot <-
  sliced_sample_tbl |>
  filter(sample_2 %in% (samples_to_inspect |> pull(sample_id))) |>
  filter(cell_number > 1) |>
  pull(file_name)

## choose layout automatically
n <- length(files_to_plot)
nc <- ceiling(sqrt(n))
nr <- ceiling(n / nc)

save_to <- "~/git_control/HPCell/dev/cellnexus-2024-scripts/samples_ncell_greater_than_one_reannotate_transformation_label_distribution.png"

# ---- open device BEFORE plotting ----
png(
  filename = save_to,
  width  = 400 * nc,   # scale with number of panels
  height = 400 * nr,
  res = 120,
  type = "cairo"
)

par(mfrow = c(nr, nc), mar = c(3,3,2,1))

purrr::walk(files_to_plot, \(f) {
  readH5AD(f, reader = "R", use_hdf5 = TRUE) |>
    assay("X") |>
    as.numeric() |>
    hist(
      main = basename(f),
      ylim = c(0, 1e5)
    )
}, .progress = TRUE)

dev.off()

# TEST THREE SAMPLEs HAVE NEGATIVE COUNTS AND PASS SCETRANSFORM. 
# Samples are used to create cellNexus paper transformation plots
sliced_sample_tbl |> filter(sample_2 %in% c("734b3e20708fa099e9af65400d3d30ba")) |> pull(file_name) |> 
  readH5AD(.x, reader = "R", use_hdf5 = T) |> assay("X") |> as.numeric() |> hist()

sliced_sample_tbl |> filter(sample_2 %in% c("734b3e20708fa099e9af65400d3d30ba")) |> 
  select(method_to_apply)

# Raw sample counts max val =0, completely negative values.
# These samples are from one dataset
sliced_sample_tbl |>
  filter(sample_2 %in% (samples_to_inspect |> pull(sample_id))) |>
  filter(cell_number > 1) |> distinct(dataset_id)
# A tibble: 1 × 1
# dataset_id                          
# <chr>                               
# 1 3f32121d-126b-4e8d-9f69-d86502d2a1b1

# Inspect the dataset counts distribution
x = readH5AD("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/h5ad/2024-07-01/3f32121d-126b-4e8d-9f69-d86502d2a1b1.h5ad",
         reader = "R", use_hdf5 = T)

x |> assay("X") |> as.numeric() |> summary()


###### Rerun failed samples ######
### Rerun all failed samples in HPCell ###
### HPCell "thorough" will rerun targets on and after the udpate one ###
### For example, input_file_function[199] <- "exp" will make HPCell rerun 199, and input_files after index 199 ###
### To Do: Issue: CREATE_PR_ISSUE ###
sample_target_tbl = tar_meta(starts_with("sct_"), store = "/vast/scratch/users/shen.m/cellNexus_target_store_2024_Jul") |> 
  filter(!is.na(error))  |> 
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

sample_target_tbl = sample_target_tbl |> bind_rows()
sample_target_tbl |> left_join( tar_meta(starts_with("sct_"), 
                      store = "/vast/scratch/users/shen.m/cellNexus_target_store_2024_Jul") |> filter(!is.na(error)) |> 
               select(name,error),
             by = c("target_name" = "name"),
             copy=T)  |> group_by(error) |> slice_head(n = 1)
# tar_workspace(sct_matrix_6aea3a48c1ea7024, store = "/vast/scratch/users/shen.m/cellNexus_target_store_2024_Jul" )
# debugonce(non_batch_variation_removal)
# non_batch_variation_removal(sce_transformed, empty_tbl, alive_tbl, doublet_tbl, 
#                             NULL, factors_to_regress = c("subsets_Mito_percent",
#                                                          "subsets_Ribo_percent"),
#                             external_path = "~/scratch/cache_temp", container_type = data_container_type)

samples_to_rerun <- sliced_sample_tbl |> filter(sample_2 %in% (sample_target_tbl |> pull(sample_id)))
samples_to_rerun |> write_parquet("~/cellxgene_curated/metadata_cellxgene_mengyuan/failed_sct_samples_to_rerun_Jul_2024.parquet")
x <-
  samples_to_rerun |> 
  pull(file_name) |> 
  str_replace("/home/users/allstaff/shen.m/scratch", "/vast/scratch/users/shen.m") |> 
  set_names(samples_to_rerun |> pull(sample_2))
f  = samples_to_rerun |> pull(method_to_apply)
tr = samples_to_rerun |> pull(feature_thresh)

job::job({
  
  library(HPCell)
  
  x |>
    initialise_hpc(
      store = "~/scratch/samples_failed_sct_to_rerun_target_store",
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
            #memory_gigabytes_required = c(75, 80, 90, 100, 120, 150), 
            memory_gigabytes_required = c(100, 110,  120, 140, 150), 
            cpus_per_task = c(2, 2, 5, 10, 20), 
            time_minutes = c(60*24, 60*24, 60*24, 60*24, 60*24,60*24),
            verbose = T
          )
        )
        
      ),
      verbosity = "summary",
      update = "thorough", 
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
    
    # Pseudobulk
    calculate_pseudobulk(target_input = "sce_transformed",
                         group_by = "cell_type_unified_ensemble") |>
    
    # metacell
    cluster_metacell(target_input = "sce_transformed",  group_by = "cell_type_unified_ensemble") |>
    
    print()
})

tar_script({
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(targets)
  library(tarchetypes)
  library(crew)
  library(crew.cluster)
  tar_option_set(
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue", 
    #debug = "annotation_tbl_light", 
    cue = tar_cue(mode = "never"), 
    controller = crew_controller_group(
      list(crew_controller_slurm(
        name = "tier_2",
        script_lines = "#SBATCH --mem 10G",
        slurm_cpus_per_task = 1,
        workers = 200,
        tasks_max = 10,
        verbose = T,
        #launch_max = 5, 
        seconds_idle = 30,
        slurm_time_minutes = 480
      ) )
    ), 
    trust_object_timestamps = TRUE
  )
  
  lighten_annotation = function(target_name, my_store ){
    annotation_tbl = tar_read_raw( target_name,  store = my_store )
    if(annotation_tbl |> is.null()) { 
      warning("this annotation is null -> ", target_name)
      return(NULL) 
    }
    
    annotation_tbl |> 
      unnest(blueprint_scores_fine) |> 
      select(.cell, blueprint_first.labels.fine, monaco_first.labels.fine, any_of("azimuth_predicted.celltype.l2"), monaco_scores_fine, contains("macro"), contains("CD4") ) |> 
      unnest(monaco_scores_fine) |> 
      select(.cell, blueprint_first.labels.fine, monaco_first.labels.fine, any_of("azimuth_predicted.celltype.l2"), contains("macro") , contains("CD4"), contains("helper"), contains("Th")) |> 
      rename(cell_ = .cell)
  }
  
  list(
    
    # The input DO NOT DELETE
    tar_target(my_store, "/vast/scratch/users/shen.m/samples_failed_sct_to_rerun_target_store", deployment = "main"),
    
    tar_target(
      target_name,
      tar_meta(
        starts_with("annotation_tbl_"), 
        store = my_store) |> 
        filter(type=="branch") |> 
        pull(name),
      deployment = "main"
    )    ,
    
    tar_target(
      annotation_tbl_light,
      lighten_annotation(target_name, my_store),
      packages = c("dplyr", "tidyr"),
      pattern = map(target_name)
    )
  )
  
  
}, script = "/vast/scratch/users/shen.m/samples_failed_sct_lighten_annotation_tbl_target.R", ask = FALSE)

job::job({
  
  tar_make(
    script = "/vast/scratch/users/shen.m/samples_failed_sct_lighten_annotation_tbl_target.R", 
    store = "/vast/scratch/users/shen.m/samples_failed_sct_lighten_annotation_tbl_target", 
    reporter = "summary"
  )
  
})


samples_to_rerun <- read_parquet("~/cellxgene_curated/metadata_cellxgene_mengyuan/failed_sct_samples_to_rerun_Jul_2024.parquet")
sample_ids <- samples_to_rerun |> pull(sample_2)
library(duckdb)
# Write annotation light
cell_metadata <- 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata.parquet')")
  )  |> 
  filter(sample_id %in% sample_ids) |>
  mutate(cell_ = paste0(cell_, "___", dataset_id)) |> 
  select(cell_, observation_joinid, contains("cell_type"), dataset_id,  self_reported_ethnicity, tissue, donor_id,  sample_id, is_primary_data, assay)


cell_annotation = 
  tar_read(annotation_tbl_light, store = "/vast/scratch/users/shen.m/samples_failed_sct_lighten_annotation_tbl_target") |> 
  dplyr::rename(
    blueprint_first_labels_fine = blueprint_first.labels.fine, 
    monaco_first_labels_fine = monaco_first.labels.fine, 
    azimuth_predicted_celltype_l2 = azimuth_predicted.celltype.l2
  ) 

cell_annotation = cell_annotation |> mutate(
  blueprint_first_labels_fine = ifelse(is.na(blueprint_first_labels_fine), "Other", blueprint_first_labels_fine),
  monaco_first_labels_fine = ifelse(is.na(monaco_first_labels_fine), "Other", monaco_first_labels_fine),
  azimuth_predicted_celltype_l2=ifelse(is.na(azimuth_predicted_celltype_l2), "Other", azimuth_predicted_celltype_l2))

cell_annotation |> arrow::write_parquet("/vast/scratch/users/shen.m/cellNexus_run/samples_failed_sct_annotation_tbl_light.parquet")

empty_tbl = tar_read(empty_tbl, store = "/vast/scratch/users/shen.m/samples_failed_sct_to_rerun_target_store") |> bind_rows() |> 
  dplyr::rename(cell_ = .cell)

alive_tbl = 
  tar_read(alive_tbl, store = "/vast/scratch/users/shen.m/samples_failed_sct_to_rerun_target_store") |>
  bind_rows() |>
  dplyr::rename(cell_ = .cell)

doublet_tbl =
  tar_read(doublet_tbl, store = "/vast/scratch/users/shen.m/samples_failed_sct_to_rerun_target_store") |>
  bind_rows() |>
  dplyr::rename(cell_ = .cell)

annotation_tbl = tar_read(annotation_tbl, store = "/vast/scratch/users/shen.m/samples_failed_sct_to_rerun_target_store") |>
  bind_rows() |>
  unnest(blueprint_scores_fine) |> 
  select(.cell, blueprint_first.labels.fine, monaco_first.labels.fine, any_of("azimuth_predicted.celltype.l2"), monaco_scores_fine, contains("macro"), contains("CD4") ) |>
  unnest(monaco_scores_fine) |> 
  select(.cell, blueprint_first.labels.fine, monaco_first.labels.fine, any_of("azimuth_predicted.celltype.l2"), contains("macro") , contains("CD4"), contains("helper"), contains("Th")) |>
  dplyr::rename(cell_ = .cell)


cell_type_concensus_tbl = tar_read(cell_type_concensus_tbl, store = "/vast/scratch/users/shen.m/samples_failed_sct_to_rerun_target_store") |>  
  bind_rows() |> 
  dplyr::rename(cell_ = .cell)

cell_type_concensus_tbl = cell_type_concensus_tbl |> mutate(cell_type_unified_ensemble = 
                                                              ifelse(is.na(cell_type_unified_ensemble),
                                                                     "Unknown",
                                                                     cell_type_unified_ensemble)) 

metacell = 
  tar_read(metacell_tbl, store = "/vast/scratch/users/shen.m/samples_failed_sct_to_rerun_target_store") |> 
  bind_rows() |> 
  dplyr::rename(cell_ = cell) |> 
  dplyr::rename_with(
    ~ stringr::str_replace(.x, "^gamma", "metacell_"),
    starts_with("gamma")
  )


cell_metadata_joined = cell_metadata |> 
  left_join(empty_tbl, copy=TRUE) |>  
  left_join(cell_type_concensus_tbl, copy=TRUE) |>
  left_join(alive_tbl, copy=TRUE) |> 
  left_join(doublet_tbl, copy=TRUE) |>
  left_join(metacell, copy=TRUE)

cell_metadata_joined = cell_metadata_joined |> as_tibble() |> 
  # Match to how pseudobulk annotations get parsed in HPCell/R/functions preprocessing_output()
  mutate(cell_type_unified_ensemble = ifelse(cell_type_unified_ensemble |> is.na(), "Unknown", cell_type_unified_ensemble)) |>
  mutate(data_driven_ensemble = ifelse(data_driven_ensemble |> is.na(), "Unknown", data_driven_ensemble))   |>
  mutate(blueprint_first_labels_fine = ifelse(blueprint_first_labels_fine |> is.na(), "Other", blueprint_first_labels_fine)) |> 
  mutate(monaco_first_labels_fine = ifelse(monaco_first_labels_fine |> is.na(), "Other", monaco_first_labels_fine)) |> 
  mutate(azimuth_predicted_celltype_l2 = ifelse(azimuth_predicted_celltype_l2 |> is.na(), "Other", azimuth_predicted_celltype_l2)) |> 
  mutate(azimuth = ifelse(azimuth |> is.na(), "Other", azimuth)) |> 
  mutate(blueprint = ifelse(blueprint |> is.na(), "Other", blueprint)) |> 
  mutate(monaco = ifelse(monaco |> is.na(), "Other", monaco))

cell_metadata_remain_untouch = tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_annotation_2024_Jul.parquet')")
) |> filter(!sample_id%in%sample_ids) |> collect()

gc()

cell_metadata_joined = cell_metadata_joined |> collect()
# Merge
cell_metadata_joined_updated = bind_rows(cell_metadata_remain_untouch,cell_metadata_joined )
gc()
cell_metadata_joined_updated |> arrow::write_parquet("/vast/scratch/users/shen.m/cellNexus_run/cell_annotation_2024_Jul_updated.parquet", compression = "zstd")




