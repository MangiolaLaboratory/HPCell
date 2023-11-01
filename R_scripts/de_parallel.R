library(tidyverse)
library(jascap)
library(targets)
library(glue)


# Get input

my_data = 
  tidySummarizedExperiment::se |> 
  tidybulk::keep_abundant()
my_data_df = 
  tibble(name = "my_data", se = list(my_data)) 
my_store = tempfile(tmpdir = ".")



my_data_df |> 
  hpcell_test_differential_abundance(store = my_store)


se =
  tidySummarizedExperiment::se |>
  tidybulk::keep_abundant()
  
  se |>
    hpcell_test_differential_abundance(~ dex + (1 | cell))