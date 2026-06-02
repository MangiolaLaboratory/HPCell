library(testthat)
library(HPCell)

# ---------------------------------------------------------------------------
# HPCell pipeline integration tests
# ---------------------------------------------------------------------------
# These tests verify that the grammar API builds HPCell objects correctly and
# that a minimal end-to-end pipeline executes without error.
# Heavy-dependency tests are guarded with skip_if_not_installed() and
# skip_on_cran() so they are skipped in environments where optional packages
# are absent.
# ---------------------------------------------------------------------------

test_that("initialise_hpc returns an HPCell object with correct structure", {
  skip_if_not_installed("Seurat")
  skip_on_cran()

  store <- file.path(tempdir(), paste0("hpcell_struct_", Sys.getpid()))
  file  <- file.path(tempdir(), paste0("test_input_", Sys.getpid(), ".rds"))

  counts <- matrix(
    rpois(300 * 80, lambda = 2),
    nrow = 300, ncol = 80,
    dimnames = list(paste0("Gene", seq_len(300)), paste0("Cell", seq_len(80)))
  )
  Seurat::CreateSeuratObject(counts = counts) |> saveRDS(file)

  hpc <- c(sample1 = file) |>
    initialise_hpc(
      gene_nomenclature  = "symbol",
      data_container_type = "seurat_rds",
      store              = store
    )

  expect_s3_class(hpc, "HPCell")
  expect_equal(hpc$initialisation$data_container_type, "seurat_rds")
  expect_false(is.null(hpc$initialisation$store))
})

test_that("HPCell grammar functions chain without error", {
  skip_if_not_installed("Seurat")
  skip_on_cran()

  store <- file.path(tempdir(), paste0("hpcell_chain_", Sys.getpid()))
  file  <- file.path(tempdir(), paste0("test_chain_", Sys.getpid(), ".rds"))

  counts <- matrix(
    rpois(300 * 80, lambda = 2),
    nrow = 300, ncol = 80,
    dimnames = list(paste0("Gene", seq_len(300)), paste0("Cell", seq_len(80)))
  )
  Seurat::CreateSeuratObject(counts = counts) |> saveRDS(file)

  hpc <- c(sample1 = file) |>
    initialise_hpc(
      gene_nomenclature  = "symbol",
      data_container_type = "seurat_rds",
      store              = store
    ) |>
    remove_empty_DropletUtils() |>
    annotate_cell_type()        |>
    remove_dead_scuttle()       |>
    score_cell_cycle_seurat()   |>
    remove_doublets_scDblFinder() |>
    calculate_pseudobulk()

  expect_s3_class(hpc, "HPCell")
})

test_that("Full HPCell pipeline runs on small Seurat dataset", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("scRNAseq")
  skip_if_not_installed("DropletUtils")
  skip_if_not_installed("Azimuth")
  skip_on_cran()

  store <- file.path(tempdir(), paste0("hpcell_full_", Sys.getpid()))
  file  <- file.path(tempdir(), paste0("test_full_", Sys.getpid(), ".rds"))

  scRNAseq::HeOrganAtlasData(tissue = "Blood", ensembl = FALSE, location = FALSE) |>
    Seurat::as.Seurat(data = NULL) |>
    saveRDS(file)

  result <- c(sample1 = file) |>
    initialise_hpc(
      gene_nomenclature  = "symbol",
      data_container_type = "seurat_rds",
      store              = store
    ) |>
    remove_empty_DropletUtils() |>
    annotate_cell_type()        |>
    remove_dead_scuttle()       |>
    score_cell_cycle_seurat()   |>
    remove_doublets_scDblFinder() |>
    calculate_pseudobulk()      |>
    evaluate_hpc()

  expect_true(!is.null(result))
})
