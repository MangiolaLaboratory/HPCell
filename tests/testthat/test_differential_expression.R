library(testthat)
library(HPCell)

se =
  tidySummarizedExperiment::se |>
  tidybulk::keep_abundant(factor_of_interest = dex) |>
  test_differential_abundance_hpc(
    ~ dex + (1 | cell)
  )

test_that("simple de", {
  se |> expect_s4_class("SummarizedExperiment")
})
