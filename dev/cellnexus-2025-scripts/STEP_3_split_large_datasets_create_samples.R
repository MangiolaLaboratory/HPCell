library(zellkonverter)
library(tibble)
library(dplyr)

# for dataset containing more than 15000 cells, create sample_chunk, the first 15000 being 1, so on so force.
# new_cell_id should be dataset_id___index___sample_chunk.
# sample_id should be dataset_id___sample_chunk.
# for dataset containing less than 15000 cells, sample_id should be dataset_id

# # Parameters
# chunk_size <- 15000
# set.seed(123)   # for reproducibility if any sampling/shuffling is introduced later
# 
# # One dataset from new census
# id <- "149b2c3f-ee11-47a7-984b-923570280bd7"
# id <- "7f7faf6b-f11d-4f07-bc1c-188a4472748d"
# path <- file.path("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/h5ad/", date)
# 
# # Read
# sce <- readH5AD(
#   file.path(path, paste0(id, ".h5ad")),
#   reader = "R",
#   use_hdf5 = TRUE
# )
# 
# n_cells <- ncol(sce)
# 
# # ---- Assign chunks ----
# # This doesnt make sense because one dataset could have more than assay, donor_id
# # Deterministic: 1:15000 = chunk 1, 15001:30000 = chunk 2, etc.
# if (n_cells > chunk_size) {
#   sample_chunk <- ceiling(seq_len(n_cells) / chunk_size)
#   
#   sample_id <- paste(id, sample_chunk, sep = "___")
#   
#   new_cell_id <- paste(
#     id,
#     seq_len(n_cells),
#     sample_chunk,
#     sep = "___"
#   )
# } else {
#   sample_chunk <- NA_integer_
#   
#   sample_id <- rep(id, n_cells)
#   
#   new_cell_id <- paste(
#     id,
#     seq_len(n_cells),
#     sep = "___"
#   )
# }
# 
# # ---- Construct dataframe ----
# df <- tibble(
#   cell_id     = colnames(sce),
#   observation_joinid = colData(sce)$observation_joinid,
#   dataset_id  = id,
#   cell_chunk       = seq_len(n_cells),
#   sample_chunk = sample_chunk,
#   sample_id   = sample_id,
#   new_cell_id = new_cell_id
# )
# 
# df


# R Script for Processing Sample Metadata in Cellxgene Data
# This script processes sample metadata related to Cellxgene datasets, focusing on Homo sapiens data.
# It filters datasets based on certain criteria like primary data, accepted assays, and large sample size thresholds.
# Additionally, it modifies cell identifiers and merges this information with related datasets to generate final outputs for further analysis.
# The script employs several R packages like arrow, targets, glue, dplyr, and more for data manipulation and storage operations.


library(arrow)
library(targets)
library(glue)
library(dplyr)
library(cellxgene.census)
library(stringr)
library(purrr)
library(duckdb)

# Change metadata_cellxgenedp_Jan_2026 when iterate to new version
result_directory = "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026"
Date <- "2025-11-08"

cellxgene_dataset = 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/dataset.parquet')")
  ) |> 
  # Why they are there?
  select(-contains("X_pca"), -contains("X_umap")) |>
  distinct()
# |> 
#   filter(dataset_id %in% c("d7f5d8d0-6150-48d7-b094-c34286ad11a1", "72955cdb-bd92-4135-aa52-21f33f9640db")) 


sample_to_cell =
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_ids_for_metadata.parquet')")
  ) 

sample_to_cell = sample_to_cell |> left_join(cellxgene_dataset, by = c("dataset_id", "donor_id", "assay"),
                                                             copy = T) |> 
  
  # Discard mismatch samples
  filter(collection_id |> is.na() |> not())

# Maybe better not to produce sample_2 further
# sample_to_cell_primary |> dplyr::count(sample_) |> pull(n) |> summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0     1.0     1.0    11.6     1.0 86999.0 

# sample_to_cell_primary_human <- sample_to_cell_primary |> 
#   select(-contains("cell_type"), -is_primary_data) |> 
#   left_join(sample_meta |> select(-contains("X_pca"), -contains("X_umap"), -is_primary_data ) |> 
#               filter(organism == "Homo sapiens") |> distinct(), 
#             by = c("sample_","dataset_id", "donor_id", "sample_heuristic"),
#             copy = T) |> 
#   select(observation_joinid, cell_, sample_, 
#          dataset_id, donor_id, sample_heuristic,
#          organism, collection_id)

# accepted_assays from census 
# accepted_assays <- read.csv("~/git_control/cellNexus/dev/census_accepted_assays_2024-07-01.csv", header=TRUE) this file was used to publish cellNexus, then Census updated assay acceptance 
# url <- "https://raw.githubusercontent.com/chanzuckerberg/cellxgene-census/d44bebd3e112ea41d00aa9b2509e2a606402c07d/docs/census_accepted_assays.csv"
# download.file(url, destfile = "./dev/census_accepted_assays_2025-01-30.csv")
# accepted_assays  <- read.csv("~/git_control/cellNexus/dev/census_accepted_assays_2025-01-30.csv")
# colnames(accepted_assays) <- c("id", "assay")
# 
# #sample_to_cell_primary_human_accepted_assay <- sample_to_cell_primary_human |> filter(assay %in% accepted_assays$assay)
# 
# large_samples <- sample_to_cell_primary |> 
#   dplyr::count(sample_, collection_id, dataset_id) |> 
#   mutate(above_threshold = n > 15000)
# 
# large_samples_collection_id <- large_samples |> ungroup() |> 
#   dplyr::count(collection_id) |> arrange(desc(n))
# 
# # function to discard nucleotide in cell_ ---------------------------------
# # cell pattern repeated across samples. 
# # Decision: use modified_cell and sample_ to split data
# 
# # drop cell ID if cell ID is a series of numbers
# # ACGT more than 5, drops
# # drop cellID if does not have special cahracter : - _
# remove_nucleotides_and_separators <- function(x) {
#   # convert integer cell ID or contain numerics surrounded by special characters to NA
#   x[str_detect(x, "^[0-9:_\\-*]+$")] <- NA
#   
#   # drop sequence having a consistent stretch of 5 characters from ACGT 
#   modified <- str_replace_all(x, "[ACGT]{5,}", "")
#   
#   #remove nucleotides surrounded by optional separators
#   modified <- str_replace_all(modified, "[:_-]{2,}", "_")
# }
# 
# # List of collection IDs
# collection_ids <- large_samples_collection_id |> collect() |> pull(collection_id) |> unique()
# 
# process_collection <- function(id) {
#   filtered_data <- sample_to_cell_primary |>
#     filter(collection_id == id) |>
#     select(cell_, sample_)
#   
#   filtered_data <- filtered_data |>
#     mutate(
#       cell_modified = sql("
#       regexp_replace(
#         regexp_replace(
#           CASE 
#             WHEN regexp_matches(cell_, '^[0-9:_\\-*]+$') 
#                  THEN NULL 
#             ELSE cell_ 
#           END,
#           '[ACGT]{5,}', ''
#         ),
#         '[:_-]{2,}', '_'
#       )
#     ")
#     )
#   
#   filtered_data
# }
# # process_collection(collection_ids[[1]])
# final_result <- map(collection_ids, process_collection, .progress = T)
# final_result <- purrr::reduce(final_result, union_all)
# 
# # conditional generating sample_2 based on whether number of cells > 15K.
# sample_to_cell_primary <- sample_to_cell_primary |> 
#   left_join(large_samples, by = c("sample_","collection_id","dataset_id"))
# 
# new_sample_df <- 
#   sample_to_cell_primary |> 
#   left_join(final_result, by = c("cell_","sample_")) |>
#   # manual adjust 
#   mutate(
#     cell_modified = ifelse(dataset_id == "b2dda353-0c96-42df-8dcd-1ea7429a6feb" & sample_ == "5951a81f1d40153bab5d2b808e384f39",
#                            "s14",
#                            cell_modified),
#     cell_modified = ifelse(dataset_id == "b2dda353-0c96-42df-8dcd-1ea7429a6feb" & sample_ == "7313173de022921da50c34ea2f87c7af",
#                            "s3",
#                            cell_modified)
#   ) |>
#   mutate(sample_2 = if_else(above_threshold,
#                             paste(sample_, cell_modified, sep = "___"),
#                             sample_)
#   )
# save result
#sample_to_cell_primary_human_accepted_assay_sample_2 |> arrow::write_parquet("~/scratch/Census_rerun/sample_to_cell_primary_human_accepted_assay_sample_2_modify.parquet")

# Load Census census_version = "2025-01-30"
# Get new datasets from ~/git_control/cellNexus/dev/STEP_1_census_dataset_ids_download_path.R
# samples <- tbl(
#   dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
#   sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/census_new_datasets_2025-01-30.parquet')")
# )
# 
# census_samples_to_download <- samples |> 
#   filter(assay %in% accepted_assays$assay) |> 
#   select(observation_joinid, donor_id,dataset_id, is_primary_data, tissue, development_stage, assay, sex, self_reported_ethnicity, disease, cell_type) |>
#   
#   left_join(sample_to_cell_primary_human_accepted_assay_sample_2 |>
#               select(observation_joinid, cell_, sample_, dataset_id,
#                      sample_heuristic, organism, collection_id, n, above_threshold, cell_modified, sample_2 ),
#             by = c("observation_joinid", "dataset_id"),
#             #relationship = "many-to-many",
#             copy = T) |>
#   
#   # remove space in the sample_2, as sample_2 will be regarded as filename 
#   mutate(sample_2 = if_else(str_detect(sample_2, " "), str_replace_all(sample_2, " ",""), sample_2)) |>
#   # Remove possible duplicates
#   distinct() 

# This is important: please make sure observation_joinid and cell_ is unique per sample (sample_2) in census_samples_to_download
sample_to_cell |> dplyr::count(observation_joinid, sample_) |> dplyr::count(n)
sample_to_cell |> dplyr::count(cell_, sample_) |> dplyr::count(n)
# If above step return not unique, go to diagnostic script: ~/git_control/cellNexus/dev/optional_diagnose_duplicate_cell_ids_in_step3.R

# Check cell count per sample, if a sample has cell more than 1e6, possibly missed a special spliter in STEP2/sample_heuristic function
sample_to_cell |> dplyr::count(sample_id) |> arrange(desc(n))

sample_to_cell |>
  filter(is_primary_data) |>
  group_by(dataset_id, sample_id, assay)  |>
  summarise(observation_joinid = list(observation_joinid), .groups = "drop") |> 
  collect() |> 
  mutate(list_length = map_dbl(observation_joinid, length)) |>
  arrow::write_parquet(glue::glue("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/{Date}_census_samples_to_download.parquet"), compression = "zstd")


# Create cell_metadata
cell_ids_for_metadata <- tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_ids_for_metadata.parquet')")
)

# # this metadata is generated in ~/git_control/cellNexus/dev/STEP_3_split_large_datasets_create_samples.R
# cell_to_refined_sample_from_Mengyuan <- tbl(
#   dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
#   sql(glue::glue("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/{Date}_census_samples_to_download.parquet')"))
# ) |>
#   select(cell_, observation_joinid, dataset_id, sample_id) 


# Establish a connection to DuckDB in memory
job::job({
  
  con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # Create views for each of the datasets in DuckDB
  #   dbExecute(con, "
  #   CREATE VIEW cell_to_refined_sample_from_Mengyuan AS
  #   SELECT cell_, observation_joinid, dataset_id, sample_id, cell_type, cell_type_ontology_term_id
  #   FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/2025-11-08_census_samples_to_download.parquet')
  # ")
  
  dbExecute(con, "
  CREATE VIEW cell_ids_for_metadata AS
  SELECT *
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_ids_for_metadata.parquet')
")
  
  dbExecute(con, "
  CREATE VIEW sample_metadata AS
  SELECT *
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/dataset.parquet')
")
  
  dbExecute(con, "
  CREATE VIEW age_days_tbl AS
  SELECT development_stage, age_days
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/age_days.parquet')
")
  
  dbExecute(con, "
  CREATE VIEW tissue_grouped AS
  SELECT tissue, tissue_groups
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/tissue_grouped.parquet')
")
  
  # Perform optimised joins within DuckDB
  copy_query <- "
COPY (
  SELECT 
    cell_ids_for_metadata.*,
    sample_metadata.*,
    age_days_tbl.age_days,
    tissue_grouped.tissue_groups,
    'cellxgene' AS atlas_id
  
  FROM cell_ids_for_metadata
  
  -- LEFT JOIN cell_ids_for_metadata
  --   ON cell_ids_for_metadata.cell_ = cell_to_refined_sample_from_Mengyuan.cell_
  --   AND cell_ids_for_metadata.observation_joinid = cell_to_refined_sample_from_Mengyuan.observation_joinid
  --   AND cell_ids_for_metadata.dataset_id = cell_to_refined_sample_from_Mengyuan.dataset_id
  
  LEFT JOIN sample_metadata
    ON sample_metadata.dataset_id = cell_ids_for_metadata.dataset_id
    AND sample_metadata.donor_id = cell_ids_for_metadata.donor_id
    AND sample_metadata.assay = cell_ids_for_metadata.assay
    
  LEFT JOIN age_days_tbl
    ON age_days_tbl.development_stage = cell_ids_for_metadata.development_stage

  LEFT JOIN tissue_grouped
    ON tissue_grouped.tissue = cell_ids_for_metadata.tissue
  
  WHERE cell_ids_for_metadata.is_primary_data = TRUE
    
) TO '/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_metadata.parquet'
(FORMAT PARQUET, COMPRESSION 'gzip');
"
  
  # Execute the final query to write the result to a Parquet file
  dbExecute(con, copy_query)
  
  # Disconnect from the database
  dbDisconnect(con, shutdown = TRUE)
  
})

# system("~/bin/rclone copy /vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_metadata.parquet box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/")


cell_metadata = tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_metadata.parquet')")
) 

