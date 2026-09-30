ensemble_annotation <- function(celltype_matrix, method_weights = NULL, override_celltype = c(), celltype_tree = NULL) {
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

add_celltype_level <- function(.data, id_col, level = 0, celltype_tree = NULL) {
  if (is.null(celltype_tree)) {
    .data_internal(immune_graph)
  }
  stopifnot(is(celltype_tree, "igraph"))
  stopifnot(igraph::is_directed(celltype_tree))

  ig_diameter = igraph::diameter(celltype_tree)
  if (level > ig_diameter) {
    stop(sprintf("The specified level (%d) exceeds the depth of the celltype tree (%d)", level, ig_diameter))
  }

  # check column exists
  id_col_str = rlang::as_string(rlang::ensym(id_col))
  if (!id_col_str %in% colnames(.data)) {
    stop(sprintf("Column '%s' not found in .data", rlang::as_string(rlang::ensym(id_col))))
  }

  # generate map
  ct_map = igraph::ego(celltype_tree, mode = "in", order = ig_diameter) |>
    sapply(\(x) {
      x = rev(x$name)
      x[min(length(x), level + 1)]
    }) |>
    setNames(igraph::V(celltype_tree)$name)

  # retain types of the matching level only
  d = igraph::distances(celltype_tree, mode = "in")
  d[is.infinite(d)] = NA
  ct_level = apply(d, 1, max, na.rm = TRUE)
  is_child = igraph::degree(celltype_tree, mode = "out") == 0
  ct_map[ct_level[ct_map] != level & !is_child] = NA_character_
  map_df = data.frame(ctypes, ct_map[ctypes])
  colnames(map_df) = c(id_col_str, sprintf("%s_L%d", id_col_str, level))

  # join and return
  .data |>
    dplyr::left_join(map_df, copy = TRUE)
}

ensemble_annotation2 <- function(celltype_matrix,celltype_tree = NULL,
                                 local_conservative_node_index=numeric(0),
                                 global_conservative_node_list=character(0),
                                 reference_overide = character(0),
                                 non_consensus_label = "not hematopoietic",
                                 return_uncertainty=FALSE){
  # prepare immune graph
  if (is.null(celltype_tree)) {
    celltype_tree  = data(HCAO_graph)
  }
  node_names = igraph::V(celltype_tree)$name
  path_to_root = lapply(node_names,function(x){
    igraph::as_ids(igraph::ego(celltype_tree,nodes = x,order = -1,mode = "in")[[1]])
  })
  names(path_to_root) = node_names
  
  # check celltype_matrix
  if (ncol(celltype_matrix) <= 1) {
    # no ensemble required
    return(celltype_matrix[,1,drop=TRUE])
  } else {
    celltype_matrix = as.matrix(celltype_matrix)
    invalid_types = setdiff(celltype_matrix, c(node_names, NA))
    if (length(invalid_types) > 0) {
      warning(sprintf("the following cell types in 'celltype_matrix' are not in the graph and will be set to NA:\n"), utils::capture.output(utils::str(invalid_types)))
    }
    celltype_matrix[celltype_matrix %in% invalid_types] = NA
  }
  
  # get consensus label
  consensus_label = do.call(rbind,apply(celltype_matrix,1,function(x){
    # check for NAs in the label
    valid_labels = x[!is.na(x)]
    num_valid_labels = length(valid_labels)
    if(num_valid_labels<=0) return(data.frame(label=NA,uncertainty=1))
    if(num_valid_labels<=1) return(data.frame(label=valid_labels[1],uncertainty=1))
    
    # generate local conservative node list
    local_conservative_node_list = x[local_conservative_node_index]
    local_conservative_node_list = local_conservative_node_list[!is.na(local_conservative_node_list)]
    
    # find lowest common ancestor
    common_ancestors = Reduce(intersect,path_to_root[valid_labels])
    common_ancestor_length = length(common_ancestors)
    if(common_ancestor_length<=0) return(data.frame(label=non_consensus_label,uncertainty=1))
    lca = common_ancestors[1]
    
    # check if path to LCA overlapped
    paths_to_lca = lapply(path_to_root[valid_labels],function(xx){
      c(setdiff(xx,common_ancestors),lca)
    })
    subset_ca = do.call(c,lapply(1:(num_valid_labels-1),function(subset_size){
      these_subset = combn(valid_labels,m = subset_size,simplify = FALSE)
      these_ca = sapply(these_subset,function(nn){
        Reduce(intersect,paths_to_lca[nn])[1]
        })
      these_path_to_lca = lapply(path_to_root[these_ca],function(xx){
        c(setdiff(xx,common_ancestors),lca)
      })
      these_path_cover_other_vertex = mapply(FUN = function(xx,vv){
        all(setdiff(valid_labels,vv)%in%xx)
      },these_path_to_lca,these_subset,SIMPLIFY = TRUE)
      if(any(these_path_cover_other_vertex)) return(these_ca[which(these_path_cover_other_vertex)[1]])
      return(character(0))
    }))
    
    # calculate uncertainty score
    lca_path_length = sapply(paths_to_lca,length)-1
    ave_lca_path_length = length(lca_path_length)/sum(1/(lca_path_length))
    uncertainty_score = ave_lca_path_length/(ave_lca_path_length+common_ancestor_length)
    
    # return consensus label based on option arguments
    if(sum(valid_labels%in%reference_overide)==1) return(data.frame(label=valid_labels[valid_labels%in%reference_overide],uncertainty=1))
    if(lca%in%c(global_conservative_node_list,local_conservative_node_list)) return(data.frame(label=lca,uncertainty=uncertainty_score))
    if(length(subset_ca)<=0) subset_ca = lca
    return(data.frame(label=subset_ca[1],uncertainty=uncertainty_score))
  },simplify = FALSE))
  
  #
  if(return_uncertainty) return(consensus_label)
  return(consensus_label$label)
}

check_improvement = function(original_label,new_label,celltype_tree = NULL){
  if (is.null(celltype_tree)) {
    celltype_tree  = data(immune_graph)
  }
  node_names = igraph::V(celltype_tree)$name
  path_to_root = lapply(node_names,function(x){
    igraph::as_ids(igraph::ego(celltype_tree,nodes = x,order = -1,mode = "in")[[1]])
  })
  names(path_to_root) = node_names
  improvement_status = mapply(FUN = function(ori.lab,new.lab){
    this_distance=NA
    if(ori.lab==new.lab){
      this_distance = 0
    }else{
      if(new.lab%in%path_to_root[[ori.lab]]){
        this_distance = which(path_to_root[[ori.lab]]==new.lab)-1
      }else{
        if(ori.lab%in%path_to_root[[new.lab]]){
          this_distance = 1-which(path_to_root[[new.lab]]==ori.lab)
        }
      } 
    }
    return(this_distance)
  },original_label,new_label)
}
  
  
  
  