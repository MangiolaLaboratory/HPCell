# This scripts test cellNexus API with new data generated in ~/git_control/cellNexus/dev/STEP_7_unify_cell_metadata.R
library(dplyr)
library(cellNexus)
library(stringr)

cache = "/vast/scratch/users/shen.m/cellNexus/"

x = get_metadata(cache_directory = cache, cloud_metadata = NULL, local_metadata = "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/metadata.v2025.1.0.0.parquet")
sce = x |>
  dplyr::filter(
    self_reported_ethnicity == "Brazilian" &
      assay |> stringr::str_like("%10x%")
  ) |> 
  get_single_cell_experiment(cache_directory = "/vast/scratch/users/shen.m/cellNexus", repository = NULL,
                             assay = c("counts", "cpm", "rank"))

# Check the number of cells per dataset
x |> dplyr::count(dataset_id)

# Check the number of cells per sample_id
x |> dplyr::count(sample_id) |> dplyr::count(n>10)

# Check whether cell_id strategt is implemented
x |> select(cell_id) |> arrange(cell_id)

# Check QC metrics
x |> dplyr::count(empty_droplet, alive, scDblFinder.class)

# Check empty droplet ratio
x |> dplyr::count(empty_droplet) |>
  collect() |> mutate(n_cells = sum(n), pct = n / n_cells * 100)

# Check alive ratio
x |> filter(!empty_droplet) |> dplyr::count(alive) |> 
  collect() |> mutate(n_cells = sum(n), pct = n / n_cells * 100)

# Check doublet ratio
x |> filter(!empty_droplet, alive) |> dplyr::count(scDblFinder.class) |> 
  collect() |> mutate(n_cells = sum(n), pct = n / n_cells * 100)
