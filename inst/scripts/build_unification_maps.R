library(tidyverse)

# load mappings between predictions and our dictionary
map_files = list.files("inst/extdata", pattern = "immune_map.+.csv", full.names = TRUE)
names(map_files) = gsub("immune_map_(.+).csv", "\\1", basename(map_files))
celltype_unification_maps = map_files |>
  lapply(read.csv)

nonimmune_cellxgene = celltype_unification_maps$cellxgene |>
  filter(!is_immune) |>
  pull("to") |>
  unique()

# harmonise to a common nomenclature
celltype_unification_maps$azimuth = celltype_unification_maps$azimuth |>
  select(from, to) |>
  dplyr::rename(
    azimuth_predicted_celltype_l2 = from,
    azimuth = to
  )
celltype_unification_maps$blueprint = celltype_unification_maps$blueprint |>
  select(from, to) |>
  dplyr::rename(
    blueprint_first_labels_fine = from,
    blueprint = to
  )
celltype_unification_maps$monaco = celltype_unification_maps$monaco |>
  select(from, to) |>
  dplyr::rename(
    monaco_first_labels_fine = from,
    monaco = to
  )
celltype_unification_maps$cellxgene = celltype_unification_maps$cellxgene |>
  select(from, to) |>
  dplyr::rename(
    cell_type = from,
    cell_type_unified = to
  )

usethis::use_data(celltype_unification_maps)
usethis::use_data(nonimmune_cellxgene)
