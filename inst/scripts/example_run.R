data(celltype_unification_maps)
data(nonimmune_cellxgene)

# unify cell types
cell_metadata = cell_metadata |>
  left_join(celltype_unification_maps$azimuth, copy = TRUE) |>
  left_join(celltype_unification_maps$blueprint, copy = TRUE) |>
  left_join(celltype_unification_maps$monaco, copy = TRUE) |>
  left_join(celltype_unification_maps$cellxgene, copy = TRUE) |>
  mutate(ensemble_joinid = paste(azimuth, blueprint, monaco, cell_type_unified, sep = "_"))

# produce the ensemble map
df_map = cell_metadata |>
  count(ensemble_joinid, azimuth, blueprint, monaco, cell_type_unified, name = "NCells") |>
  as_tibble() |>
  mutate(
    cellxgene = if_else(cell_type_unified %in% nonimmune_cellxgene, "non immune", cell_type_unified),
    data_driven_ensemble = ensemble_annotation(cbind(azimuth, blueprint, monaco), override_celltype = c("non immune", "nkt", "mast")),
    cell_type_unified_ensemble = ensemble_annotation(cbind(azimuth, blueprint, monaco, cellxgene), method_weights = c(1, 1, 1, 2), override_celltype = c("non immune", "nkt", "mast")),
    cell_type_unified_ensemble = case_when(
      cell_type_unified_ensemble == "non immune" & cellxgene == "non immune" ~ cell_type_unified,
      cell_type_unified_ensemble == "non immune" & cellxgene != "non immune" ~ "other",
      .default = cell_type_unified_ensemble
    ),
    is_immune = !cell_type_unified_ensemble %in% nonimmune_cellxgene
  ) |>
  select(
    ensemble_joinid,
    data_driven_ensemble,
    cell_type_unified_ensemble,
    is_immune
  )

# use map to perform cell type ensemble
cell_metadata = cell_metadata |>
  left_join(df_map, by = join_by(ensemble_joinid), copy = TRUE)