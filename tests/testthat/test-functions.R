se = 
  tidySummarizedExperiment::se |> 
  tidybulk::keep_abundant()

test_that("de analysis single", {
  
  se |> 
    hpcell_test_differential_abundance(~ dex + (1 | cell))
})

test_that("de analysis multi", {
  

  tibble(name = c("my_data_1", "my_data_2"), data = list(se, se)) |> 
      hpcell_map_test_differential_abundance(~ dex + (1 | cell), data)
  

})


test_that("split", {
  
  tibble(name = c("my_data_1"), data = list(se)) |> 
    mutate(split = 5) |> 
    map_split_se_by_gene(data, split) |> 
    nrow() |> 
    expect_equal(5)

})