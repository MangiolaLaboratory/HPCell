# Description:
# This script clean up and generate the ultimate metadata to ship to cellNexus. It reads metadata and dataset-specific information,
# cleans and renames columns, and writes updated data back to disk. The process involves
# connecting to databases in memory, executing SQL queries, and handling data in both
# Parquet and HDF5 formats.

library(duckdb)
library(dbplyr)
library(dplyr)
library(tidyr)
library(data.table)
library(HDF5Array)
library(SummarizedExperiment)
library(tidySingleCellExperiment)
library(stringr)
library(targets)
library(purrr)
library(arrow)


# Add low confidence ethnicity and imputed ethnicity labels to metadata. Both data are from Ning via email
lowConf_ethnicity_df <- zellkonverter::readH5AD("/vast/projects/cellxgene_curated/cellNexus/sce_relabel.h5ad", reader = "R", use_hdf5 = T) |>
  colData() |> as_tibble() |>
  mutate(low_confidence_ethnicity = ifelse(ethnicity_relabel == "LowConfidenceLabel", TRUE, FALSE) |> as.character()) |>
  select(sample_id, ethnicity_flagging_score = score, low_confidence_ethnicity = low_confidence_ethnicity)

imputed_ethnicity_df <- zellkonverter::readH5AD("/vast/projects/cellxgene_curated/cellNexus/adata_unlabelled_with_predictions.h5ad", reader = "R", use_hdf5 = T)|>
  colData() |> as_tibble() |>
  select(sample_id, imputed_ethnicity = ethnicity_predictions) |> 
  mutate(imputed_ethnicity = as.character(imputed_ethnicity))

# lowConf_ethnicity_df |> arrow::write_parquet("/vast/projects/cellxgene_curated/cellNexus/lowConf_ethnicity_df.parquet")
# imputed_ethnicity_df |> arrow::write_parquet("/vast/projects/cellxgene_curated/cellNexus/imputed_ethnicity_df.parquet")

job::job({
  
  duckdb_write_parquet <- function(.tbl_sql, path, con) {
    
    sql_tbl <- 
      .tbl_sql |>
      sql_render()
    
    # zstd 15 compresses faster than brotli for binary/scientific datasets, whereas brotli reduce could save 100Mb
    sql_call <- glue::glue("COPY ({sql_tbl}) TO '{path}' (FORMAT PARQUET, COMPRESSION 'brotli')")
    
    res <- dbExecute(con, sql_call)
    
    return(res)
  }
  
  # Single DuckDB connection: do the heavy transforms in SQL (avoid read/write/read on 50M+ rows)
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  raw_path <- "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_5_1_filtered_missing_cells_mengyuan.parquet"  # MODIFY HERE: Metadata input parquet path
  
  DBI::dbExecute(con, glue::glue("
  CREATE VIEW cell_metadata_raw AS
  SELECT *
  FROM read_parquet({DBI::dbQuoteString(con, raw_path)}, union_by_name=true);
  "))
  
  raw_cols <- DBI::dbGetQuery(con, "SELECT * FROM cell_metadata_raw LIMIT 0") |> names()
  
  explicit_drop <- c()
  # explicit_drop <- c(
  #   "cell_",
  #   "cell__1",
  #   "dataset_id_1",
  #   "dataset_id_1_1",
  #   "cell__2",
  #   "cell__3",
  #   "dataset_id_2",
  #   "dataset_id_3",
  #   "sample_id_1",
  #   "sample_id_2",
  #   "sample_placeholder",
  #   "cell_type_unified_ensemble_1",
  #   "cell_type_1",
  #   "dataset_id_2",
  #   "observation_joinid_1",
  #   "self_reported_ethnicity_1",
  #   "donor_id_1",
  #   "assay_1",
  #   "blueprint_first_labels_fine_1",
  #   "azimuth_predicted_celltype_l2_1",
  #   "monaco_first_labels_fine_1",
  #   "dataset_id_3",
  #   "atlas_id_1",
  #   "tissue_1",
  #   "is_primary_data_1",
  #   "cell_type_ontology_term_id_1",
  #   "azimuth",
  #   "blueprint",
  #   "monaco",
  #   "alive_1",
  #   "cell_id_1",
  #   "dataset_id_4",
  #   "X_umap1",
  #   "X_umap2",
  #   "observation_originalid",
  #   "subsets_Mito_sum",
  #   "subsets_Mito_detected",
  #   "file_id_cellNexus_single_cell_1",
  #   "ensemble_joinid",
  #   "cell_type_unified",
  #   "data_driven_ensemble"
  # )
  
  pattern_drop <- c(
    grep("^scores", raw_cols, value = TRUE),
    grep("coarse$", raw_cols, value = TRUE)
  )
  
  drop_cols <- intersect(unique(c(explicit_drop, pattern_drop)), raw_cols)
  
  int_cast_cols <- intersect(
    unique(
      c(
        "feature_count",
        "nFeature_expressed_in_sample",
        "cell_count",
        grep("metacell_", raw_cols, value = TRUE),
        grep("_chunk", raw_cols, value = TRUE),
        grep("subsets_", raw_cols, value = TRUE)
      )
    ),
    raw_cols
  )
  
  chr_cast_cols <- intersect(c("published_at", "revised_at"), raw_cols)
  
  sql_id <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
  
  # Remove originals that we re-add under new names
  base_keep <- setdiff(
    raw_cols,
    c(
      drop_cols,
      "alive",
      "atlas_id",
      "blueprint_first_labels_fine",
      "monaco_first_labels_fine",
      "azimuth_predicted_celltype_l2",
      "cell_id",
      "new_cell_id"
    )
  )
  
  select_exprs <- purrr::map_chr(base_keep, function(col) {
    col_id <- sql_id(col)
    if (col %in% int_cast_cols) {
      glue::glue("CAST({col_id} AS INTEGER) AS {col_id}")
    } else if (col %in% chr_cast_cols) {
      glue::glue("CAST({col_id} AS VARCHAR) AS {col_id}")
    } else {
      col_id
    }
  })
  
  # alive: NA -> FALSE
  if ("alive" %in% raw_cols) {
    select_exprs <- c(select_exprs, glue::glue("COALESCE({sql_id('alive')}, FALSE) AS {sql_id('alive')}"))
  } else {
    select_exprs <- c(select_exprs, glue::glue("FALSE AS {sql_id('alive')}"))
  }
  
  # Rename annotation columns
  select_exprs <- c(
    select_exprs,
    if ("blueprint_first_labels_fine" %in% raw_cols) {
      glue::glue("{sql_id('blueprint_first_labels_fine')} AS {sql_id('cell_annotation_blueprint_singler')}")
    } else {
      glue::glue("NULL::VARCHAR AS {sql_id('cell_annotation_blueprint_singler')}")
    },
    if ("monaco_first_labels_fine" %in% raw_cols) {
      glue::glue("{sql_id('monaco_first_labels_fine')} AS {sql_id('cell_annotation_monaco_singler')}")
    } else {
      glue::glue("NULL::VARCHAR AS {sql_id('cell_annotation_monaco_singler')}")
    },
    if ("azimuth_predicted_celltype_l2" %in% raw_cols) {
      glue::glue("{sql_id('azimuth_predicted_celltype_l2')} AS {sql_id('cell_annotation_azimuth_l2')}")
    } else {
      glue::glue("NULL::VARCHAR AS {sql_id('cell_annotation_azimuth_l2')}")
    }
  )
  
  # new_cell_id -> cell_id as first column (drop original cell_id entirely)
  cell_id_expr <- if ("new_cell_id" %in% raw_cols) {
    glue::glue("{sql_id('new_cell_id')} AS {sql_id('cell_id')}")
  } else {
    glue::glue("NULL::VARCHAR AS {sql_id('cell_id')}")
  }
  select_exprs <- c(cell_id_expr, select_exprs)
  
  select_sql <- paste(select_exprs, collapse = ",\n    ")
  
  DBI::dbExecute(con, glue::glue("
  CREATE OR REPLACE VIEW cell_metadata AS
  SELECT
    {DBI::SQL(select_sql)}
  FROM cell_metadata_raw
  WHERE dataset_id NOT IN ('99950e99-2758-41d2-b2c9-643edcdf6d82', '9fcb0b73-c734-40a5-be9c-ace7eea401c9');
  "))
  
  DBI::dbExecute(con, "
  CREATE OR REPLACE VIEW sample_celltype_count AS
  SELECT
    sample_id,
    cell_type_unified_ensemble,
    CAST(COUNT(*) AS INTEGER) AS \".aggregated_cells\"
  FROM cell_metadata
  WHERE empty_droplet = FALSE
    AND alive = TRUE
    AND \"scDblFinder.class\" != 'doublet'
  GROUP BY sample_id, cell_type_unified_ensemble;
  ")
  
  gc()
  
  dbExecute(con, "
  CREATE VIEW lowConf_ethnicity_df AS
  SELECT 
    *
  FROM read_parquet('/vast/projects/cellxgene_curated/cellNexus/lowConf_ethnicity_df.parquet')
")
  
  dbExecute(con, "
  CREATE VIEW imputed_ethnicity_df AS
  SELECT 
    *
  FROM read_parquet('/vast/projects/cellxgene_curated/cellNexus/imputed_ethnicity_df.parquet')
")
  
  # Perform left join and save to parquet
  # MODIFY HERE: output metadata parquet path and atlas_id
  copy_query <- "
  COPY (
    SELECT
      cell_metadata.*,
      lowConf_ethnicity_df.ethnicity_flagging_score,
      lowConf_ethnicity_df.low_confidence_ethnicity,
      sample_celltype_count.\".aggregated_cells\",
      COALESCE(imputed_ethnicity_df.imputed_ethnicity, cell_metadata.self_reported_ethnicity) AS imputed_ethnicity, -- Use imputed_ethnicity if present
      'cellxgene_2024/0.2.1' AS atlas_id
      
    FROM cell_metadata
    
    LEFT JOIN lowConf_ethnicity_df
      ON cell_metadata.sample_id = lowConf_ethnicity_df.sample_id
    
    LEFT JOIN imputed_ethnicity_df
      ON cell_metadata.sample_id = imputed_ethnicity_df.sample_id
    
    LEFT JOIN sample_celltype_count
      ON cell_metadata.sample_id = sample_celltype_count.sample_id AND cell_metadata.cell_type_unified_ensemble = sample_celltype_count.cell_type_unified_ensemble
      
    
    

  ) TO '/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/metadata.2.2.1.parquet'
  (FORMAT PARQUET, COMPRESSION 'zstd');
  "
  
  # Execute the final query to write the result to a Parquet file
  dbExecute(con, copy_query)
  
  # Disconnect from the database
  dbDisconnect(con, shutdown = TRUE)
  
  print("Done.")
  
  
})

x = tbl(dbConnect(duckdb::duckdb(), dbdir = ":memory:"),  
        sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/metadata.2.2.1.parquet')") ) # MODIFY HERE: input metadata parquet path

# Split cell_metadata to cellnexus_metadata, original census_metadata, and metacell_metadata (host Rshiny on smaller file)
# ---- Split: read metadata.x.y.z.parquet once, write smaller derivative Parquets ----
# This avoids loading the full metadata into R memory and avoids duplicate-column issues in Parquet writes.
job::job({
  
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  
  input_metadata <- "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/metadata.2.2.1.parquet" # MODIFY HERE: input metadata parquet path
  out_dir <- "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan"
  
  DBI::dbExecute(
    con,
    glue::glue(
      "
    CREATE OR REPLACE VIEW metadata AS
    SELECT *
    FROM read_parquet({DBI::dbQuoteString(con, input_metadata)}, union_by_name=true);
    "
    )
  )
  
  cols <- DBI::dbGetQuery(con, "SELECT * FROM metadata LIMIT 0") |> names()
  sql_id <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
  
  # Strip sample annotation from cellnexus annotation doesn't save too much (less than 3Mb), thus keep in one.
  remove_cols <- c(
    "cell_type", "cell_type_ontology_term_id", "data_driven_ensemble", "ensemble_joinid",
    "observation_originalid", "assay", "assay_ontology_term_id", "development_stage", "development_stage_ontology_term_id",
    "disease", "disease_ontology_term_id", "donor_id", "is_primary_data", "organism", "organism_ontology_term_id",
    "self_reported_ethnicity", "self_reported_ethnicity_ontology_term_id",
    "sex", "sex_ontology_term_id", "tissue", "tissue_ontology_term_id", "citation",
    "collection_id", "dataset_version_id", "default_embedding", "published_at", "raw_data_location",
    "revised_at", "primary_cell_count", "schema_version", "tissue_type", "title",
    "tombstone", "x_approximate_distribution", "explorer_url", "cell_count", "feature_count", 
    "filesize", "filetype", "mean_genes_per_cell", "suspension_type", "url"
  )
  
  # CellNexus metadata (smaller file for Shiny): drop heavy / internal columns by name patterns
  drop_cellnexus <- unique(c(
    intersect(remove_cols, cols),
    cols[grepl("metacell", cols)]
  ))
  keep_cellnexus <- setdiff(cols, drop_cellnexus)
  select_cellnexus <- paste(sql_id(keep_cellnexus), collapse = ", ")
  
  # MODIFY HERE: output cellnexus metadata parquet path
  DBI::dbExecute(
    con,
    glue::glue(
      "
      COPY (
        SELECT {DBI::SQL(select_cellnexus)}
        FROM metadata
      )
      TO {DBI::dbQuoteString(con, file.path(out_dir, 'cellnexus_metadata.2.2.1.parquet'))}
      (FORMAT PARQUET, COMPRESSION 'brotli');
      "
    )
  )
  
  # Original census-like metadata subset (stable columns)
  census_cols <- intersect(
    c(
      "observation_joinid", "dataset_id", "sample_id", "cell_type",
      "cell_type_ontology_term_id", "assay", "assay_ontology_term_id", "development_stage", "development_stage_ontology_term_id",
      "disease", "disease_ontology_term_id", "donor_id", "is_primary_data", "organism", "organism_ontology_term_id",
      "self_reported_ethnicity", "self_reported_ethnicity_ontology_term_id",
      "sex", "sex_ontology_term_id", "tissue", "tissue_ontology_term_id",
      "data_driven_ensemble", "ensemble_joinid", "observation_originalid",  "citation",
      "collection_id", "dataset_version_id", "default_embedding", "published_at", "raw_data_location",
      "revised_at", "primary_cell_count", "schema_version", "tissue_type", "title",
      "tombstone", "x_approximate_distribution", "explorer_url", "cell_count", "feature_count", 
      "filesize", "filetype", "mean_genes_per_cell", "suspension_type", "url"
    ),
    cols
  )
  select_census <- paste(sql_id(census_cols), collapse = ", ")
  
  # MODIFY HERE: output census metadata parquet path
  DBI::dbExecute(
    con,
    glue::glue(
      "
    COPY (
      SELECT {DBI::SQL(select_census)}
      FROM metadata
    )
    TO {DBI::dbQuoteString(con, file.path(out_dir, 'census_cell_metadata.2.2.1.parquet'))}
    (FORMAT PARQUET, COMPRESSION 'brotli');
    "
    )
  )
  
  # # Metacell metadata subset
  # metacell_cols <- unique(c("cell_id", "sample_id", "dataset_id", cols[grepl("metacell", cols)]))
  # metacell_cols <- intersect(metacell_cols, cols)
  # select_metacell <- paste(sql_id(metacell_cols), collapse = ", ")
  # 
  # # MODIFY HERE: output metacell metadata parquet path
  # DBI::dbExecute(
  #   con,
  #   glue::glue(
  #     "
  #   COPY (
  #     SELECT {DBI::SQL(select_metacell)}
  #     FROM metadata
  #   )
  #   TO {DBI::dbQuoteString(con, file.path(out_dir, 'metacell_metadata.2.2.1.parquet'))}
  #   (FORMAT PARQUET, COMPRESSION 'brotli');
  #   "
  #   )
  # )
  
  print("Done.")
})


# (Optional) Check whether cellnexus_metadata parquet can be optimised further
# source("~/git_control/cellNexus/dev/data_optimisation_script.R")

