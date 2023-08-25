test_that("de analysis", {
  
  tidySummarizedExperiment::se |> 
    tidybulk::keep_abundant() |> 
    hpcell_test_differential_abundance()
})
