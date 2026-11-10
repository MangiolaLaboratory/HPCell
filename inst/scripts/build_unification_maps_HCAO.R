library(tidyverse)

# load mappings between predictions and our dictionary
map_files = list.files("inst/extdata/", pattern = "immune_map.+._to_CL.csv", full.names = TRUE)
names(map_files) = gsub("immune_map_(.+)_to_CL.csv", "\\1", basename(map_files))
HCAO_celltype_unification_maps = map_files |>
  lapply(read.csv)

# harmonise to the CL nomenclature
HCAO_celltype_unification_maps$azimuth = HCAO_celltype_unification_maps$azimuth |>
  select(from, to) |>
  dplyr::rename(
    azimuth_predicted_celltype_l2 = from,
    azimuth = to
  )
HCAO_celltype_unification_maps$blueprint = HCAO_celltype_unification_maps$blueprint |>
  select(from, to) |>
  dplyr::rename(
    blueprint_first_labels_fine = from,
    blueprint = to
  )
HCAO_celltype_unification_maps$monaco = HCAO_celltype_unification_maps$monaco |>
  select(from, to) |>
  dplyr::rename(
    monaco_first_labels_fine = from,
    monaco = to
  )
HCAO_celltype_unification_maps$cellxgene = HCAO_celltype_unification_maps$cellxgene |>
  select(from, to) |>
  dplyr::rename(
    cell_type_ontology_term_id = from,
    cellxgene = to
  )

usethis::use_data(HCAO_celltype_unification_maps)

