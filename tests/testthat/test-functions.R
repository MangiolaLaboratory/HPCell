test_that("de analysis single", {
  
  tidySummarizedExperiment::se |> 
    tidybulk::keep_abundant() |> 
    hpcell_test_differential_abundance(~ dex + (1 | cell))
})

test_that("de analysis multi", {
  
  se = 
    tidySummarizedExperiment::se |> 
    tidybulk::keep_abundant()
  
  tibble(name = c("my_data_1", "my_data_2"), data = list(se, se)) |> 
      hpcell_map_test_differential_abundance(~ dex + (1 | cell), data)
  

})
