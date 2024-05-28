test_that("get_census_data() downloads the correct data type in a defined directory", {
  path <- tempfile()

  census <- open_soma(census_version = "2023-12-15")
  metadata <-
    census$get("census_data")$get("homo_sapiens")$get("obs")
  selected_columns <-
    c(
      "assay",
      "disease",
      "donor_id",
      "sex",
      "self_reported_ethnicity",
      "tissue",
      "development_stage"
    )
  sample <-
    metadata$read(column_names = selected_columns)$concat() |>
    as.data.frame() |>
    distinct() |>
    head(1)

  # downloads data
  sample |> get_census_data(
    save_path = path,
    experiment = "sce",
    data_output_type = "rds"
  )
  filename <-
    list.files(path, pattern = "\\.rds$", full.names = TRUE)
  file <- readRDS(filename)

  # downloads data in a correct type
  file |> expect_s4_class("SingleCellExperiment")

  file |> rownames() |>
    length() |>
    expect_gt(1)
})

test_that("rename_features() renames ensembl features to gene names", {
  path <- tempfile()

  census <- open_soma(census_version = "2023-12-15")
  metadata <-
    census$get("census_data")$get("homo_sapiens")$get("obs")
  selected_columns <-
    c(
      "assay",
      "disease",
      "donor_id",
      "sex",
      "self_reported_ethnicity",
      "tissue",
      "development_stage"
    )
  sample <-
    metadata$read(column_names = selected_columns)$concat() |> as.data.frame() |> distinct() |> head(1)

  sample |> get_census_data(
    save_path = path,
    experiment = "sce",
    data_output_type = "rds"
  )
  filename <-
    list.files(path, pattern = "\\.rds$", full.names = TRUE)

  file <- readRDS(filename)

  file |> rownames() |>
    stringr::str_like("^ENSG*") |>
    all() |>
    expect_false()
})

test_that("HPCell pipeline works for the Seurat test data", {
  tempdir <- tempdir()
  file <- file.path(tempdir, "file.rds")
  scRNAseq::HeOrganAtlasData(tissue = "Blood",
                             ensembl = FALSE,
                             location = FALSE) |>
    as.Seurat(data = NULL) |>
    saveRDS(file)
  
  # Running the pipeline
  run_targets_pipeline(
    input_data = file,
    tissue = "pbmc",
    sample_column = "Tissue",
    filter_empty_droplets = TRUE,
    store = tempdir(),
    data_container_type = "seurat_rds"
  )
  
  preprocessing_output <-
    tar_read(preprocessing_output_S, store = tempdir())[[1]]
  
  preprocessing_output[[]] |> colnames() |>
    stringr::str_which(pattern = "subsets_Mito_percent|subsets_Ribo_percent") |>
    expect_length(2)
})
