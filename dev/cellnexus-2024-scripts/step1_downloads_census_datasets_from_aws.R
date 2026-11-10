# This script identifies CELLxGENE Census datasets stable version from 2024-07-01, and downloads
# their corresponding .h5ad files from AWS S3.
#
# Workflow:
# 1. For every dataset, generate AWS S3 download commands pointing to:
#       /vast/scratch/users/shen.m/test_cellnexus_reproducibility/h5ad/<release_date>/
# 2. Write all download commands into a text file.
# 3. Create a GNU Parallel bash script that downloads all new datasets efficiently.
library(cellxgene.census)
library(dplyr)

# Retrieve previous census stable version dataset_id 
census <- open_soma(census_version = "2024-07-01")
# Identify organisms available in the Census
census_data <- census$get("census_data")
org_name <- names(census_data$members)[grepl("homo", names(census_data$members), ignore.case = TRUE)]

metadata <- census_data$get(org_name)$get("obs")

selected_columns <- c('assay', 'assay_ontology_term_id','disease', 'disease_ontology_term_id',
                      'donor_id', 'sex', 'sex_ontology_term_id', 'self_reported_ethnicity', 'self_reported_ethnicity_ontology_term_id',
                      'tissue', 'tissue_ontology_term_id', 'tissue_type', 'development_stage', 'development_stage_ontology_term_id',
                      'is_primary_data','dataset_id','observation_joinid', 'suspension_type',
                      "cell_type", "cell_type_ontology_term_id")

samples <- metadata$read(column_names = selected_columns,
                         value_filter = "is_primary_data == 'TRUE'")$concat()

# Get new datasets
samples <- samples |> as.data.frame() |> distinct()  |>
  # Add organism
  mutate(organism = org_name)

# saved <- samples |> arrow::write_parquet(glue::glue("/vast/scratch/users/shen.m/test_cellnexus_reproducibility/census_new_datasets_{date}.parquet"),
#                                          compression = "zstd")
# 

# Set the base path where files will be downloaded
date <- "2024-07-01"
path = file.path("/vast/scratch/users/shen.m/test_cellnexus_reproducibility/h5ad/", date)
h5ads.uri = get_census_version_directory() |> tibble::rownames_to_column("version") |> filter(version == date) |> pull(h5ads.uri)

if (!dir.exists(path)) dir.create(path, recursive = T)

# Generate the dataset download commands
# Only new datasets will be downloaded to path
dataset_ids_path <- samples |> distinct(dataset_id) |> 
  #head(2 ) |>
  mutate(download_path =  paste0("aws s3 cp --no-sign-request ",  h5ads.uri,  dataset_id, ".h5ad ", path)) |> select(download_path)

# Path where the script will be saved
output_file_path <- glue::glue(
  "/vast/scratch/users/shen.m/test_cellnexus_reproducibility/h5ad/{date}_census_dataset_ids_download_path.txt"
)

write.table(
  dataset_ids_path,  
  file = output_file_path,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Define the content of the bash script
bash_script <- glue::glue("

#!/bin/bash

module load awscli
module load parallel

# Disable all /dev/tty output from GNU Parallel
export PARALLEL_NO_TTY=1

# Config to boost download speed
aws configure set default.s3.max_bandwidth 70MB/s

# Max 100 parallel HTTPs connections for one file
# Limit AWS CLI internal concurrency PER download (prevents connection storms)
aws configure set default.s3.max_concurrent_requests 8

# Path to the file containing AWS S3 copy commands
COMMAND_FILE='{output_file_path}'

# Number of parallel downloads. Try lower the variable here to avoid TCP limits.
# At most 1-2 simultaneous downloads (prevents firewall congestion)
PARALLEL_DOWNLOADS=2

# Execute the download commands in parallel
cat $COMMAND_FILE | parallel -j $PARALLEL_DOWNLOADS --eta --bar --plain --no-notice
")

# Write the bash script to a file
writeLines(bash_script, "/vast/scratch/users/shen.m/test_cellnexus_reproducibility/parallel_download.sh")

# Change file permission to make it executable
system("chmod +x /vast/scratch/users/shen.m/test_cellnexus_reproducibility/parallel_download.sh")

# Execute the bash script
system("/vast/scratch/users/shen.m/test_cellnexus_reproducibility/parallel_download.sh")
