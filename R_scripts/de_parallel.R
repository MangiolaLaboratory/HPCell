library(tidyverse)
library(HPCell)
library(targets)
library(glue)
library(tictoc)

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
    hpcell_test_differential_abundance(
      ~ dex + (1 | cell),
      computing_resources = crew.cluster::crew_controller_slurm(
        name = "slurm",
        slurm_memory_gigabytes_per_cpu = 5,
        slurm_cpus_per_task = 2,
        workers = 200,
        verbose = T
      )
      )
