library(testthat)
library(HPCell)

# ---------------------------------------------------------------------------
# Differential expression tests
# ---------------------------------------------------------------------------

test_that("test_differential_abundance_hpc returns a SummarizedExperiment", {
  skip_if_not_installed("tidySummarizedExperiment")
  skip_if_not_installed("tidybulk")
  skip_if_not_installed("lme4")
  skip_on_cran()

  se <- tidySummarizedExperiment::se |>
    tidybulk::keep_abundant(factor_of_interest = dex) |>
    test_differential_abundance_hpc(~ dex + (1 | cell))

  expect_s4_class(se, "SummarizedExperiment")
})
