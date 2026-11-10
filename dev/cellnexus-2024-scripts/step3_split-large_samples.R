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
result_directory = "/vast/scratch/users/shen.m/test_cellnexus_reproducibility"
# # Sample metadata
# sample_meta <- tar_read(metadata_dataset_id_common_sample_columns, store = glue("{result_directory}/_targets"))
# sample_meta |> arrow::write_parquet("~/scratch/Census/sample_meta.parquet", compression = "zstd")

# Sample to cell link
# sample_to_cell <- tar_read(metadata_dataset_id_cell_to_sample_mapping, store = glue("{result_directory}/_targets"))
# sample_to_cell_primary <- sample_to_cell |> filter(is_primary_data == TRUE)
# sample_to_cell_primary |> arrow::write_parquet("~/scratch/Census/sample_to_cell_primary.parquet", compression = "zstd")

sample_meta = 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/sample_metadata.parquet')")
  )

sample_to_cell_primary = 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_ids_for_metadata.parquet')")
  )

sample_to_cell_primary_human <- sample_to_cell_primary |> 
  left_join(sample_meta |> filter(organism == "Homo sapiens"), 
            by = c("sample_","dataset_id", "donor_id", "is_primary_data", "sample_heuristic"),
            copy = T) |> 
  select(observation_joinid, cell_, sample_, 
         dataset_id, donor_id, is_primary_data, sample_heuristic,
         organism, tissue, development_stage, assay, collection_id,
         sex, self_reported_ethnicity, disease)
gc()

# accepted_assays from census 
# accepted_assays <- read.csv("~/git_control/cellNexus/dev/census_accepted_assays_2024-07-01.csv", header=TRUE) this file was used to publish cellNexus, then Census updated assay acceptance 
# url <- "https://raw.githubusercontent.com/chanzuckerberg/cellxgene-census/d44bebd3e112ea41d00aa9b2509e2a606402c07d/docs/census_accepted_assays.csv"
# download.file(url, destfile = "./dev/census_accepted_assays_2025-01-30.csv")
accepted_assays  <- read.csv("~/git_control/cellNexus/dev/census_accepted_assays_2024-07-01.csv")
colnames(accepted_assays) <- c("id", "assay")

sample_to_cell_primary_human_accepted_assay <- sample_to_cell_primary_human |> filter(assay %in% accepted_assays$assay)

large_samples <- sample_to_cell_primary_human_accepted_assay |> 
  dplyr::count(sample_, assay, collection_id, dataset_id) |> 
  mutate(above_threshold = n > 15000)

large_samples_collection_id <- large_samples |> ungroup() |> 
  dplyr::count(collection_id) |> arrange(desc(n))

# function to discard nucleotide in cell_ ---------------------------------
# cell pattern repeated across samples. 
# Decision: use modified_cell and sample_ to split data

# drop cell ID if cell ID is a series of numbers
# ACGT more than 5, drops
# drop cellID if does not have special cahracter : - _
remove_nucleotides_and_separators <- function(x) {
  # convert integer cell ID or contain numerics surrounded by special characters to NA
  x[str_detect(x, "^[0-9:_\\-*]+$")] <- NA
  
  # drop sequence having a consistent stretch of 5 characters from ACGT 
  modified <- str_replace_all(x, "[ACGT]{5,}", "")
  
  #remove nucleotides surrounded by optional separators
  modified <- str_replace_all(modified, "[:_-]{2,}", "_")
}

# List of collection IDs
collection_ids <- large_samples_collection_id |> collect() |> pull(collection_id) 

gc()
process_collection <- function(id) {
  filtered_data <- sample_to_cell_primary_human_accepted_assay |>
    filter(collection_id == id) |>
    collect() |>
    select(cell_, sample_)
  
  #filtered_data <- filtered_data |> mutate(cell_modified = remove_nucleotides_and_separators(cell_))
  filtered_data$cell_modified <- remove_nucleotides_and_separators(filtered_data |> pull(cell_))
  filtered_data
}

final_result <- map(collection_ids, process_collection, .progress = T)
final_result <- reduce(final_result, union_all)

# conditional generating sample_2 based on whether number of cells > 15K.
sample_to_cell_primary_human_accepted_assay <- sample_to_cell_primary_human_accepted_assay |> 
  left_join(large_samples, by = c("sample_", "assay","collection_id","dataset_id"))

sample_to_cell_primary_human_accepted_assay_sample_2 <- 
  sample_to_cell_primary_human_accepted_assay |> 
  left_join(final_result, by = c("cell_","sample_"), copy = TRUE) |>
  # manual adjust 
  mutate(
    cell_modified = ifelse(dataset_id == "b2dda353-0c96-42df-8dcd-1ea7429a6feb" & sample_ == "5951a81f1d40153bab5d2b808e384f39",
                           "s14",
                           cell_modified),
    cell_modified = ifelse(dataset_id == "b2dda353-0c96-42df-8dcd-1ea7429a6feb" & sample_ == "7313173de022921da50c34ea2f87c7af",
                           "s3",
                           cell_modified)
  ) |>
  mutate(sample_2 = if_else(above_threshold,
                            paste(sample_, cell_modified, sep = "___"),
                            sample_)
  )
# save result
write_parquet_to_parquet = function(data_tbl, output_parquet, compression = "gzip") {
  
  # Establish connection to DuckDB in-memory database
  con_write <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # Register `data_tbl` within the DuckDB connection (this doesn't load it into memory)
  duckdb::duckdb_register(con_write, "data_tbl_view", data_tbl)
  
  # Use DuckDB's COPY command to write `data_tbl` directly to Parquet with compression
  copy_query <- paste0("
  COPY data_tbl_view TO '", output_parquet, "' (FORMAT PARQUET, COMPRESSION '", compression, "');
  ")
  
  # Execute the COPY command
  dbExecute(con_write, copy_query)
  
  # Unregister the temporary view
  duckdb::duckdb_unregister(con_write, "data_tbl_view")
  
  # Disconnect from the database
  dbDisconnect(con_write, shutdown = TRUE)
}

sample_to_cell_primary_human_accepted_assay_sample_2 |> write_parquet_to_parquet("~/scratch/Census_rerun/sample_to_cell_primary_human_accepted_assay_sample_2_modify.parquet")

gc()

# Load Census census_version = "2024-07-01"
census <- open_soma(census_version = "2024-07-01")
metadata <- census$get("census_data")$get("homo_sapiens")$get("obs")
selected_columns <- c('assay', 'disease', 'donor_id', 'sex', 'self_reported_ethnicity', 'tissue', 'development_stage','is_primary_data','dataset_id','observation_joinid',
                      "cell_type", "cell_type_ontology_term_id")
samples <- metadata$read(column_names = selected_columns,
                         value_filter = "is_primary_data == 'TRUE'")$concat()
samples <- samples |> as.data.frame() |> distinct() 
samples |> write_parquet("~/scratch/Census/census_samples_701.parquet")

######## READ
#sample_to_cell_primary_human_accepted_assay_sample_2 <- arrow::read_parquet("~/scratch/Census_rerun/sample_to_cell_primary_human_accepted_assay_sample_2.parquet")
samples <- 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('~/scratch/Census/census_samples_701.parquet')")
  )

sample_to_cell_primary_human_accepted_assay_sample_2 <- tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('~/scratch/Census_rerun/sample_to_cell_primary_human_accepted_assay_sample_2_modify.parquet')")
)

census_samples_to_download <- samples |> 
  left_join(sample_to_cell_primary_human_accepted_assay_sample_2,
            by = c("observation_joinid", "dataset_id"),
            relationship = "many-to-many",
            copy=TRUE) |>
  # Use annotation from census
  select(-donor_id.y,
         -is_primary_data.y,
         -tissue.y,
         -development_stage.y,
         -assay.y,
         -sex.y,
         -self_reported_ethnicity.y,
         -disease.y) |>
  rename(donor_id = donor_id.x,
         is_primary_data = is_primary_data.x,
         assay = assay.x,
         disease = disease.x,
         sex = sex.x,
         self_reported_ethnicity = self_reported_ethnicity.x,
         tissue = tissue.x,
         development_stage = development_stage.x
  ) |> 
  #as_tibble() |>
  # remove space in the sample_2, as sample_2 will be regarded as filename 
  mutate(sample_2 = if_else(str_detect(sample_2, " "), str_replace_all(sample_2, " ",""), sample_2))

# For query purpose 
census_samples_to_download |> write_parquet_to_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/census_samples_to_download_MODIFIED.parquet")


# This is important: please make sure observation_joinid and cell_ is unique per sample (sample_2) in census_samples_to_download
census_samples_to_download |> dplyr::count(observation_joinid, sample_2) |> dplyr::count(n)
census_samples_to_download |> dplyr::count(cell_, sample_2) |> dplyr::count(n)

# light version
census_samples_to_download |> group_by(dataset_id, sample_2)  |> 
  summarise(observation_joinid = list(observation_joinid), .groups = "drop") |> as_tibble() |> mutate(list_length = map_dbl(observation_joinid, length)) |>
  arrow::write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/census_samples_to_download_groups_MODIFIED.parquet")

# Establish a connection to DuckDB in memory
job::job({
  
  con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # Create views for each of the datasets in DuckDB
  dbExecute(con, "
  CREATE VIEW cell_to_refined_sample_from_Mengyuan AS
  SELECT cell_, observation_joinid, dataset_id, sample_2 AS sample_id, cell_type, cell_type_ontology_term_id
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/census_samples_to_download_MODIFIED.parquet')
")
  
  dbExecute(con, "
  CREATE VIEW cell_ids_for_metadata AS
  SELECT cell_, observation_joinid, dataset_id, sample_, donor_id
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_ids_for_metadata.parquet')
")
  
  dbExecute(con, "
  CREATE VIEW sample_metadata AS
  SELECT *
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/sample_metadata.parquet')
")
  
  dbExecute(con, "
  CREATE VIEW age_days_tbl AS
  SELECT development_stage, age_days
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/age_days.parquet')
")
  
  dbExecute(con, "
  CREATE VIEW tissue_grouped AS
  SELECT tissue, tissue_groups
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/tissue_grouped.parquet')
")
  
  # Perform optimised joins within DuckDB
  copy_query <- "
COPY (
  SELECT 
    cell_to_refined_sample_from_Mengyuan.cell_,
    cell_to_refined_sample_from_Mengyuan.observation_joinid,
    cell_to_refined_sample_from_Mengyuan.dataset_id,
    cell_to_refined_sample_from_Mengyuan.sample_id,
    cell_to_refined_sample_from_Mengyuan.cell_type,
    cell_to_refined_sample_from_Mengyuan.cell_type_ontology_term_id,
    sample_metadata.*,
    age_days_tbl.age_days,
    tissue_grouped.tissue_groups
  
  FROM cell_to_refined_sample_from_Mengyuan
  
  LEFT JOIN cell_ids_for_metadata
    ON cell_ids_for_metadata.cell_ = cell_to_refined_sample_from_Mengyuan.cell_
    AND cell_ids_for_metadata.observation_joinid = cell_to_refined_sample_from_Mengyuan.observation_joinid
    AND cell_ids_for_metadata.dataset_id = cell_to_refined_sample_from_Mengyuan.dataset_id
    
  LEFT JOIN sample_metadata
    ON cell_ids_for_metadata.sample_ = sample_metadata.sample_
    AND cell_ids_for_metadata.donor_id = sample_metadata.donor_id
    AND cell_ids_for_metadata.dataset_id = sample_metadata.dataset_id
    
  LEFT JOIN age_days_tbl
    ON age_days_tbl.development_stage = sample_metadata.development_stage

  LEFT JOIN tissue_grouped
    ON tissue_grouped.tissue = sample_metadata.tissue
    
) TO '/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata.parquet'
(FORMAT PARQUET, COMPRESSION 'gzip');
"
  
  # Execute the final query to write the result to a Parquet file
  dbExecute(con, copy_query)
  
  # Disconnect from the database
  dbDisconnect(con, shutdown = TRUE)
  
})

# system("~/bin/rclone copy /vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata.parquet box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/")


cell_metadata = tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata.parquet')")
) 

cell_metadata |> distinct(sample_id) |>dplyr::count()
