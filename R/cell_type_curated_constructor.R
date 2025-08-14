# Define the generic function
#' @export
celltype_consensus_constructor <- function(input_hpc, 
                                         target_input = "data_object", 
                                         target_output = "cell_type_concensus_tbl", 
                                         target_annotation = "annotation_tbl",
                                         annotation_unified_names = c("azimuth", "blueprint", "monaco", "cellxgene"),
                                         celltype_unification_list = NULL,
                                         nonimmune_cellxgene = NULL,
                                         ...) {
  UseMethod("celltype_consensus_constructor")
}

#' @importFrom purrr map
#' 
#' @export
celltype_consensus_constructor.HPCell <- function(input_hpc, 
                                                target_input = "data_object", 
                                                target_output = "cell_type_concensus_tbl",
                                                target_annotation = "annotation_tbl",
                                                annotation_unified_names = c("azimuth", "blueprint", "monaco", "cellxgene"),
                                                ...) {
  
  input_hpc |> 
    
    hpc_iterate(
      target_output = target_output, 
      user_function = cell_type_ensembl_harmonised |> quote(),
      input_read_RNA_assay = target_input |> is_target(), 
      annotation_label_transfer_tbl = target_annotation |> is_target(),
      available_maps = annotation_unified_names,
      ...
    )
}

#' Harmonize Cell Types Across Datasets
#'
#' This function integrates and harmonizes cell type annotations across multiple
#' datasets by applying predefined unification maps and cell type labels.
#' It uses a combination of transferred annotations and predefined maps to
#' produce a consensus on cell type identities.
#'
#' @param input_read_RNA_assay SingleCellExperiment or Seurat object containing RNA assay data.
#' @param annotation_label_transfer_tbl A tibble with annotation label transfer data.
#' @param celltype_unification_maps A list containing mapping data frames for different sources
#'             (e.g., Azimuth, Blueprint, Monaco, and cellxgene). Default is `NULL`.
#'             If `NULL`, it retrieves default maps stored in HPCell.
#' @param nonimmune A character vector specifying non-immune cell types.
#'             Default is `NULL`. If `NULL`, it retrieves default non-immune types from HPCell.
#' @param available_maps A character vector of cell type annotation sources to include in the ensemble annotation process.
#'   Supported values include `"azimuth"`, `"blueprint"`, `"monaco"`, and `"cellxgene"`. 
#'   By default, it uses all of the annotations.
#' @return A tibble of the input SummarizedExperiment metadata enriched with unified cell type annotations
#'         and additional classification details.
#'         
#' @importFrom dplyr left_join select rename mutate count as_tibble if_else case_when
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom purrr map
#' @importFrom tidyr unnest
#' @export
cell_type_ensembl_harmonised <- function(input_read_RNA_assay,
                                         annotation_label_transfer_tbl = NULL, 
                                         celltype_unification_maps = NULL, 
                                         nonimmune = NULL,
                                         available_maps = c("azimuth", "blueprint", "monaco", "cellxgene")
                                         ) {
  
  # Handle missing input
  if (input_read_RNA_assay |> is.null()) return(NULL)
  
  # Handle empty annotation_tbl
  if (nrow(annotation_label_transfer_tbl) == 0 ) return(NULL)
  
  # Use pre-generated celltype_unification_maps (list) and 
  #   nonimmune_cellxgene (character vector) from Dharmesh
  if (is.null(celltype_unification_maps)) 
    celltype_unification_maps <- HPCell::celltype_unification_maps
  if (is.null(nonimmune)) 
    nonimmune <- HPCell::nonimmune_cellxgene
  
  # get cell_metadata
  try({
    if (inherits(annotation_label_transfer_tbl, "tbl_df")){
      input_read_RNA_assay <- input_read_RNA_assay |>
        left_join(annotation_label_transfer_tbl, by = ".cell")
    }
  }, silent = TRUE)
  
  # Get metadata 
  if (inherits(input_read_RNA_assay, "Seurat")) {
    input_read_RNA_assay <-  input_read_RNA_assay[[]] |> as.data.frame() |> 
      rownames_to_column(var = ".cell") |> 
      dplyr::rename(
        blueprint_first_labels_fine = blueprint_first.labels.fine, 
        blueprint_first_labels_coarse = blueprint_first.labels.coarse,
        monaco_first_labels_fine = monaco_first.labels.fine,
        monaco_first_labels_coarse = monaco_first.labels.coarse
      ) 
  } else if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    # Rename and unnest annotation_tbl
    input_read_RNA_assay <- input_read_RNA_assay |>  SummarizedExperiment::colData() |> as.data.frame() |> 
      rownames_to_column(var = ".cell") |> 
      dplyr::rename(
        blueprint_first_labels_fine = blueprint_first.labels.fine, 
        blueprint_first_labels_coarse = blueprint_first.labels.coarse,
        monaco_first_labels_fine = monaco_first.labels.fine,
        monaco_first_labels_coarse = monaco_first.labels.coarse
      ) 
  }
  
  # Sometimes, sce does not have azimuth annotation
  input_read_RNA_assay <- input_read_RNA_assay |> 
    mutate(azimuth_predicted_celltype_l2 = ifelse(!("azimuth_predicted.celltype.l2" %in% names(input_read_RNA_assay)),
                                                  NA, 
                                                  azimuth_predicted.celltype.l2)) |>
    unnest(blueprint_scores_fine) |> 
    select(.cell, any_of(c("observation_joinid", "observation_originalid",
           "donor_id", "dataset_id", "sample_id", "cell_type")),
           blueprint_first_labels_fine, monaco_first_labels_fine, 
           blueprint_first_labels_coarse, monaco_first_labels_coarse,
           any_of("azimuth_predicted_celltype_l2"), monaco_scores_fine, contains("macro"), contains("CD4") ) |> 
    unnest(monaco_scores_fine) |> 
    select(.cell, any_of(c("observation_joinid", "observation_originalid",
                           "donor_id", "dataset_id", "sample_id", "cell_type")),
           blueprint_first_labels_fine, monaco_first_labels_fine, 
           blueprint_first_labels_coarse, monaco_first_labels_coarse,
           any_of("azimuth_predicted_celltype_l2"), contains("macro") , contains("CD4"), contains("helper"), contains("Th"))
  
  
  # If cellxgene is available, calculate ensemble annotations accordingly
  has_cellxgene <- "cellxgene" %in% available_maps
  
  # Set method weights depending on availability of cellxgene
  if (has_cellxgene) {
    method_weights <- c(1, 1, 1, 2)
  } else {
    method_weights <- c(1, 1, 1, 0)  # 0 weight for missing cellxgene
  }
  
  # Unify cell types
  input_read_RNA_assay <- input_read_RNA_assay |>
    left_join(celltype_unification_maps$azimuth, copy = TRUE) |>
    left_join(celltype_unification_maps$blueprint,  copy = TRUE) |>
    left_join(celltype_unification_maps$monaco,  copy = TRUE) |>
    mutate(ensemble_joinid = paste(azimuth, blueprint, monaco, sep = "_"))
  
  # Conditionally join cellxgene only if it exists
  if ("cellxgene" %in% available_maps) input_read_RNA_assay <- input_read_RNA_assay |> 
    left_join(celltype_unification_maps$cellxgene, copy = TRUE) |> 
    mutate(ensemble_joinid = paste(ensemble_joinid, cell_type_unified, sep = "_"))
  
  # Produce the ensemble map
  df_map <- input_read_RNA_assay |>
    dplyr::count(across(all_of(c("azimuth", "blueprint", "monaco", "cell_type_unified", "ensemble_joinid"))), name = "NCells") |>
    as_tibble() |>
    mutate(cellxgene = if (has_cellxgene) {
      if_else(cell_type_unified %in% nonimmune_cellxgene, "non immune", 
              cell_type_unified)
    } else {NA_character_},
      data_driven_ensemble = ensemble_annotation(cbind(azimuth, blueprint, monaco), 
                                                 override_celltype = c("non immune", "nkt", "mast")),
      cell_type_unified_ensemble = ensemble_annotation(cbind(azimuth, blueprint, monaco, cellxgene), 
                                                       method_weights = method_weights, 
                                                       override_celltype = c("non immune", "nkt", "mast")),
      cell_type_unified_ensemble = if (has_cellxgene) {
        case_when(
          cell_type_unified_ensemble == "non immune" & cellxgene == "non immune" ~ cell_type_unified,
          cell_type_unified_ensemble == "non immune" & cellxgene != "non immune" ~ "other",
          TRUE ~ cell_type_unified_ensemble
        )
      } else {
        case_when(
          cell_type_unified_ensemble == "non immune" ~ "other",
          TRUE ~ cell_type_unified_ensemble
        )
      },
      is_immune = !cell_type_unified_ensemble %in% nonimmune
    ) |>
    select(
      ensemble_joinid,
      data_driven_ensemble,
      cell_type_unified_ensemble,
      is_immune
    )
  
  # Use map to perform cell type ensemble
  input_read_RNA_assay <- input_read_RNA_assay |>
    left_join(df_map, by = "ensemble_joinid", copy = TRUE)
  
  return(input_read_RNA_assay)
}


#' Ensemble Annotation for Cell Type Identification
#'
#' This function creates an ensemble annotation for cell types by utilizing a voting mechanism
#' across different methods. It leverages a hierarchy of cell types, method-specific weights,
#' and an option to override certain cell types to derive a consensus classification.
#'
#' @param celltype_matrix A matrix or data frame where columns represent different annotation
#'        methods for cell types. Each element in the matrix represents a cell type determined by
#'        each method.
#' @param method_weights Optional numeric vector or matrix specifying weights for each method.
#'        If not provided, equal weights are used. If provided as a vector, it should match the
#'        number of methods (columns of celltype_matrix).
#' @param override_celltype A character vector of cell types that should override the voting
#'        process if they appear. This can be used to set certain cell types as non-negotiable
#'        when they are detected by any method.
#' @param celltype_tree An igraph object representing the hierarchy of cell types. If NULL,
#'        a default graph named "immune_graph" from the global environment is used.
#'
#' @return A vector representing the consensus cell type for each row in the input `celltype_matrix`.
#'
#' @export
ensemble_annotation <- function(celltype_matrix, method_weights = NULL, 
                                override_celltype = c(), celltype_tree = NULL) {
  if (is.null(celltype_tree)) {
    celltype_tree <- get("immune_graph")
  }
  
  stopifnot(is(celltype_tree, "igraph"))
  stopifnot(igraph::is_directed(celltype_tree))
  stopifnot(is.matrix(celltype_matrix) | is.data.frame(celltype_matrix))
  
  node_names = igraph::V(celltype_tree)$name
  
  # check override_celltype nodes are present
  missing_nodes = setdiff(override_celltype, node_names)
  if (!is.null(missing_nodes) & length(missing_nodes) > 0) {
    missing_nodes = paste(missing_nodes, collapse = ", ")
    stop(sprintf("the following nodes in 'override_celltype' not found in 'celltype_tree': %s", utils::capture.output(utils::str(missing_nodes))))
  }
  
  # check celltype_matrix
  if (ncol(celltype_matrix) == 1) {
    # no ensemble required
    return(celltype_matrix)
  } else {
    celltype_matrix = as.matrix(celltype_matrix)
    invalid_types = setdiff(celltype_matrix, c(node_names, NA))
    if (length(invalid_types) > 0) {
      warning(sprintf("the following cell types in 'celltype_matrix' are not in the graph and will be set to NA:\n"), utils::capture.output(utils::str(invalid_types)))
    }
    celltype_matrix[celltype_matrix %in% invalid_types] = NA
  }
  
  # check method_weights
  if (is.null(method_weights)) {
    method_weights = matrix(1, ncol = ncol(celltype_matrix), nrow = nrow(celltype_matrix))
  } else if (is.vector(method_weights)) {
    if (ncol(celltype_matrix) != length(method_weights)) {
      stop("the number of columns in 'celltype_matrix' should match the length of 'method_weights'")
    }
    method_weights = matrix(rep(method_weights, each = nrow(celltype_matrix)), nrow = nrow(celltype_matrix))
  } else if (is.matrix(method_weights) | is.data.frame(method_weights)) {
    if (ncol(celltype_matrix) != ncol(method_weights)) {
      stop("the number of columns in 'celltype_matrix' and 'method_weights' should be equal")
    }
    method_weights = as.matrix(method_weights)
  }
  method_weights = method_weights / rowSums(method_weights)
  
  # create vote matrix
  vote_matrix = Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = c(nrow(celltype_matrix), length(node_names)), dimnames = list(rownames(celltype_matrix), node_names))
  for (i in seq_len(ncol(celltype_matrix))) {
    locmat = cbind(seq_len(nrow(celltype_matrix)), as.numeric(factor(celltype_matrix[, i], levels = node_names)))
    missing = is.na(locmat[, 2])
    vote_matrix[locmat[!missing, ]] = vote_matrix[locmat[!missing, ]] + method_weights[!missing, i]
  }
  
  # propagate vote to children
  d = apply(!is.infinite(igraph::distances(celltype_tree, mode = "out")), 2, as.numeric)
  d = as(d, "sparseMatrix")
  vote_matrix_children = Matrix::tcrossprod(vote_matrix, Matrix::t(d))
  
  # propagate vote to parent
  d = igraph::distances(celltype_tree, mode = "in")
  d = 1 / (2^d) - 0.1 # vote halved at each subsequent ancestor
  diag(d)[igraph::degree(celltype_tree, mode = "in") > 0 & igraph::degree(celltype_tree, mode = "out") == 0] = 0
  diag(d) = diag(d) * 0.9 # prevent leaf nodes from being selected when trying to identify upstream ancestor (works for any number in the interval (0.5, 1))
  vote_matrix_parent = Matrix::tcrossprod(vote_matrix, Matrix::t(d))
  
  # assess votes and identify common ancestors for ties
  vote_matrix_children = apply(vote_matrix_children, 1, \(x) x[x > 0], simplify = FALSE)
  vote_matrix_parent = apply(vote_matrix_parent, 1, \(x) x[x > 0], simplify = FALSE)
  ensemble = mapply(\(children, parents) {
    # override condition
    override_node = intersect(override_celltype, names(children))
    if (length(override_node) > 0) {
      return(override_node[1])
    }
    
    # maximum votes
    children = names(children)[children == max(children)]
    if (length(children) == 1) {
      return(children)
    } else {
      # lowest ancestor with the maximum votes
      parents = names(parents)[parents == max(parents)]
      if (length(parents) == 1) {
        return(parents)
      } else {
        return(NA)
      }
    }
  }, vote_matrix_children, vote_matrix_parent)
  
  return(ensemble)
}
