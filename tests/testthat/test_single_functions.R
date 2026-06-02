library(testthat)
library(HPCell)

# ---------------------------------------------------------------------------
# Unit tests for individual HPCell functions
# ---------------------------------------------------------------------------
# Tests are ordered from lightest (no heavy deps) to heaviest.
# All tests that load Seurat, Bioconductor data, or run models are guarded
# with skip_if_not_installed() and skip_on_cran().
# ---------------------------------------------------------------------------

# --- is_target -----------------------------------------------------------

test_that("is_target converts a string to a name symbol", {
  result <- is_target("my_target")
  expect_true(is.name(result))
  expect_equal(as.character(result), "my_target")
})

test_that("is_target returns NULL for NULL input", {
  expect_null(is_target(NULL))
})

test_that("is_target errors on non-character input", {
  expect_error(is_target(42))
})

# --- empty_droplet_threshold (lightweight counts matrix) -----------------

test_that("empty_droplet_threshold returns a tibble with empty_droplet column", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("EnsDb.Hsapiens.v86")
  skip_on_cran()

  set.seed(1)
  counts <- matrix(
    rpois(500 * 100, lambda = 1),
    nrow = 500, ncol = 100,
    dimnames = list(paste0("Gene", seq_len(500)), paste0("Cell", seq_len(100)))
  )
  obj <- Seurat::CreateSeuratObject(counts = counts)

  result <- HPCell:::empty_droplet_threshold(
    obj,
    feature_nomenclature  = "symbol",
    RNA_feature_threshold = 2
  )

  expect_s3_class(result, "tbl_df")
  expect_true("empty_droplet" %in% colnames(result))
  expect_true(nrow(result) > 0)
})

# --- empty_droplet_id (requires DropletUtils) ----------------------------

test_that("empty_droplet_id returns a tibble with expected columns", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("DropletUtils")
  skip_if_not_installed("EnsDb.Hsapiens.v86")
  skip_on_cran()

  set.seed(1)
  counts <- matrix(
    rpois(500 * 150, lambda = 1),
    nrow = 500, ncol = 150,
    dimnames = list(paste0("Gene", seq_len(500)), paste0("Cell", seq_len(150)))
  )
  obj <- Seurat::CreateSeuratObject(counts = counts)

  result <- HPCell:::empty_droplet_id(
    obj,
    total_RNA_count_check = -Inf,
    feature_nomenclature  = "symbol"
  )

  expect_s3_class(result, "tbl_df")
  expect_true("empty_droplet" %in% colnames(result))
})

# --- doublet_identification -----------------------------------------------

test_that("doublet_identification returns a tibble with scDblFinder.class", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("scDblFinder")
  skip_on_cran()

  set.seed(1)
  counts <- matrix(
    rpois(300 * 120, lambda = 2),
    nrow = 300, ncol = 120,
    dimnames = list(paste0("Gene", seq_len(300)), paste0("Cell", seq_len(120)))
  )
  obj <- Seurat::CreateSeuratObject(counts = counts)

  result <- HPCell:::doublet_identification(obj)

  expect_s3_class(result, "tbl_df")
  expect_true("scDblFinder.class" %in% colnames(result))
})

# --- cell_cycle_scoring ---------------------------------------------------

test_that("cell_cycle_scoring returns a tibble with Phase column", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("EnsDb.Hsapiens.v86")
  skip_on_cran()

  set.seed(1)
  counts <- matrix(
    rpois(300 * 80, lambda = 2),
    nrow = 300, ncol = 80,
    dimnames = list(paste0("Gene", seq_len(300)), paste0("Cell", seq_len(80)))
  )
  obj <- Seurat::CreateSeuratObject(counts = counts)

  result <- HPCell:::cell_cycle_scoring(
    obj,
    feature_nomenclature = "symbol"
  )

  expect_s3_class(result, "tbl_df")
  expect_true("Phase" %in% colnames(result))
  expect_true(all(result$Phase %in% c("G1", "S", "G2M")))
})

# --- annotation_label_transfer -------------------------------------------

test_that("annotation_label_transfer returns a tibble", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("celldex")
  skip_if_not_installed("SingleR")
  skip_if_not_installed("EnsDb.Hsapiens.v86")
  skip_on_cran()

  set.seed(1)
  counts <- matrix(
    rpois(500 * 80, lambda = 2),
    nrow = 500, ncol = 80,
    dimnames = list(paste0("Gene", seq_len(500)), paste0("Cell", seq_len(80)))
  )
  obj <- Seurat::CreateSeuratObject(counts = counts)

  result <- HPCell:::annotation_label_transfer(
    obj,
    feature_nomenclature = "symbol"
  )

  expect_s3_class(result, "tbl_df")
  expect_true(".cell" %in% colnames(result))
})

# --- alive_identification ------------------------------------------------

test_that("alive_identification returns a tibble with mitochondrial columns", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("scater")
  skip_if_not_installed("EnsDb.Hsapiens.v86")
  skip_on_cran()

  set.seed(1)
  counts <- matrix(
    rpois(500 * 80, lambda = 2),
    nrow = 500, ncol = 80,
    dimnames = list(paste0("Gene", seq_len(500)), paste0("Cell", seq_len(80)))
  )
  obj <- Seurat::CreateSeuratObject(counts = counts)

  result <- HPCell:::alive_identification(
    obj,
    feature_nomenclature = "symbol"
  )

  expect_s3_class(result, "tbl_df")
  expect_true("alive" %in% colnames(result) ||
              "subsets_Mito_percent" %in% colnames(result))
})
