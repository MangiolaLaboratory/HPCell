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
