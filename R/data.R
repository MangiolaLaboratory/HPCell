#' Dummy HPC Dataset
#'
#' This dataset, named `dummy_hpc`, is a synthetic example dataset used for demonstrating
#' high-performance computing (HPC) data processing and analysis techniques.
#'
#' @format A data frame with multiple rows and columns representing synthetic HPC data. Each column represents a different variable in the dataset.
#'
#' @details
#' The `dummy_hpc` dataset is created for educational and demonstration purposes. It includes
#' simulated data points that resemble typical HPC workload metrics. The dataset can be used to
#' showcase various data processing, analysis, and visualisation techniques in the context of HPC.
#'
#' @usage
#' data(dummy_hpc)
#'
#' @examples
#' # Load the dataset
#' data(dummy_hpc)
#'
#' # Display the first few rows of the dataset
#' head(dummy_hpc)
#'
#' # Example analysis: summary statistics
#' summary(dummy_hpc)
#'
#' @keywords datasets
#' @docType data
"dummy_hpc"

#' 
#' This dataset contains Ensembl gene IDs, external gene names, and chromosome names
#' retrieved using the biomaRt package.
#' 
#' @format A data frame map of ensembl_gene_id, external_gene_name and chromosome_name
#' 
#' @usage 
#' data(ensembl_genes_biomart)
#' 
#' @source biomaRt::getBM()
#' 
#' @keywords datasets
#' @docType data
"ensembl_genes_biomart"

#' CellChatDB.human database
#'
#' A curated human ligand–receptor interaction database provided by the CellChat package.
#'
#' This object is typically used as input to the CellChat pipeline. It contains signaling pathway data 
#' for cell-cell communication analysis.
#'
#' @format A list with multiple elements, each representing different parts of the signaling network.

#' @usage
#' data(CellChatDB.human)
#' 
#' @source CellChat::CellChatDB.human
#' @docType data
"CellChatDB.human"

#' Cell-type unification mapping tables
#'
#' A named list of data frames mapping cell-type labels from different
#' annotation references (Azimuth, Blueprint, Monaco, CellxGene) to a
#' unified cell-type vocabulary used by `cell_type_ensembl_harmonised()`.
#'
#' @format A named list with elements `azimuth`, `blueprint`, `monaco`, and
#'   optionally `cellxgene`. Each element is a data frame with at least two
#'   columns: the original label and the unified label.
#'
#' @usage
#' data(celltype_unification_maps)
#'
#' @keywords datasets
#' @docType data
"celltype_unification_maps"

#' Immune cell-type hierarchy graph
#'
#' An `igraph` directed graph encoding the hierarchical relationships between
#' immune cell types. Used by `ensemble_annotation()` and `add_celltype_level()`
#' to resolve consensus cell-type calls across annotation methods.
#'
#' @format An `igraph` object where each vertex is a cell type and edges
#'   represent parent–child relationships in the hierarchy.
#'
#' @usage
#' data(immune_graph)
#'
#' @keywords datasets
#' @docType data
"immune_graph"

#' Non-immune CellxGene cell-type labels
#'
#' A character vector of cell-type labels from the CellxGene corpus that are
#' classified as non-immune. Used by `cell_type_ensembl_harmonised()` to
#' distinguish immune from non-immune populations during consensus annotation.
#'
#' @format A character vector.
#'
#' @usage
#' data(nonimmune_cellxgene)
#'
#' @keywords datasets
#' @docType data
"nonimmune_cellxgene"

#' HCAO immune cell-type hierarchy graph
#'
#' An `igraph` directed graph of immune cell types from the Human Cell Atlas
#' Ontology (HCAO), with vertices named by Cell Ontology IDs (e.g. `CL_0000785`)
#' plus a root-level `"not hematopoietic"` node. Used by `ensemble_annotation2()`
#' to resolve consensus cell-type calls across annotation methods.
#'
#' @format An `igraph` object with 94 vertices and 92 directed edges, where
#'   each edge points from a parent cell type to a child cell type.
#'
#' @usage
#' data(HCAO_graph)
#'
#' @source Built by `inst/scripts/build_graph_HCAO.R` from
#'   `inst/extdata/HCAO_immune_graph_edges.csv` and
#'   `inst/extdata/HCAO_immune_graph_vertices.csv`.
#' @keywords datasets
#' @docType data
"HCAO_graph"

#' HCAO cell-type unification mapping tables
#'
#' A named list of data frames mapping cell-type labels from different
#' annotation references (Azimuth, Blueprint, Monaco, CellxGene) to Cell
#' Ontology IDs that are vertices of [HCAO_graph].
#'
#' @format A named list with elements `azimuth`, `blueprint`, `cellxgene`, and
#'   `monaco`. Each element is a data frame with two columns: the original
#'   reference label and the unified Cell Ontology ID.
#'
#' @usage
#' data(HCAO_celltype_unification_maps)
#'
#' @source Built by `inst/scripts/build_unification_maps_HCAO.R`.
#' @keywords datasets
#' @docType data
"HCAO_celltype_unification_maps"
