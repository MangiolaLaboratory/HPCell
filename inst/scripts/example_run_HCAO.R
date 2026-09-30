library(tidyverse)

data(HCAO_celltype_unification_maps)
data(HCAO_graph)

source(here::here("R/functions_consensus.R"))

# load cellNexus metadata and only keep high quality cells
### This will download ~1GB of data
cellNexus_metadata = cellNexus::get_metadata() |> 
  cellNexus::keep_quality_cells()

# load census cell_type ontology
census_metadata = cellNexus:::get_census_metadata("2024-07-01")
con = dbplyr::remote_con(cellNexus_metadata)

duckdb::duckdb_register_arrow(con, "census_metadata", census_metadata)

# join census cell type ontology label back to cellnexus
cellNexus_metadata = cellNexus_metadata |>
  dplyr::left_join(tbl(con, "census_metadata") |> 
                     dplyr::select(observation_joinid, dataset_id, 
                                   cell_type, cell_type_ontology_term_id))


# unify cell types
HCAO_celltype_unification_maps$cellxgene = 
  HCAO_celltype_unification_maps$cellxgene |>
  mutate(cell_type_ontology_term_id = gsub(pattern = "_",replacement = ":",cell_type_ontology_term_id,fixed = TRUE))

cellNexus_metadata = cellNexus_metadata |>
  left_join(HCAO_celltype_unification_maps$azimuth, by = join_by(cell_annotation_azimuth_l2==azimuth_predicted_celltype_l2), copy = TRUE) |>
  left_join(HCAO_celltype_unification_maps$blueprint, by = join_by(cell_annotation_blueprint_singler==blueprint_first_labels_fine), copy = TRUE) |>
  left_join(HCAO_celltype_unification_maps$monaco, by = join_by(cell_annotation_monaco_singler==monaco_first_labels_fine), copy = TRUE) |>
  left_join(HCAO_celltype_unification_maps$cellxgene, by = join_by(cell_type_ontology_term_id), copy=TRUE)|>
  mutate(cellxgene = case_when(
    is.na(cellxgene) ~ "not hematopoietic",
    .default = cellxgene
  ),
    ensemble_joinid = paste(azimuth, blueprint, monaco, cellxgene, sep = "__"))

df_map = cellNexus_metadata |>
  count(ensemble_joinid, azimuth, blueprint, monaco, cellxgene, name = "NCells") |>
  collect()

# monoply exception for blueprint's macrophage
macrophage_exception = c("CL_0000235","CL_0000863","CL_0000861","CL_0000890")

# no-reference exception for cellxgene's mast cell and natural killer T cell
mast_and_nkt_exception = c("CL_0000814","CL_0000097")

# get data-driven consensus and ensemble consensus
df_map = df_map |>
  mutate(
    data_driven_ensemble = ensemble_annotation2(
      cbind(azimuth, blueprint, monaco),
      celltype_tree = HCAO_graph,
      reference_overide = macrophage_exception)
  )|>
  mutate(
    cell_type_unified_ensemble_HCAO = ensemble_annotation2(
      cbind(data_driven_ensemble, cellxgene),
      local_conservative_node_index = c(1),
      celltype_tree = HCAO_graph,
      reference_overide = mast_and_nkt_exception)
  )|>
  mutate(
    cell_type_unified_ensemble_HCAO = case_when(
      cell_type_unified_ensemble_HCAO == "not hematopoietic" & cellxgene != "not hematopoietic" ~ "other",
      .default = cell_type_unified_ensemble_HCAO
    ),
    is_immune_HCAO = ! cell_type_unified_ensemble_HCAO %in% c("not hematopoietic","other")
  )

# use map to perform cell type ensemble
cellNexus_metadata = cellNexus_metadata |>
  left_join(df_map |>
              select(ensemble_joinid,
                     NCells,
                     data_driven_ensemble,
                     cell_type_unified_ensemble_HCAO,
                     is_immune_HCAO), 
            by = join_by(ensemble_joinid), copy = TRUE)
