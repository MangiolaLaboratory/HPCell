# Input: seurat, output: nested tibble of variable features
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat ScaleData
#' @importFrom Seurat RunPCA
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom Seurat VariableFeatures
#' @importFrom glue glue
seurat_to_variable_features_by_cell_type = function(counts, assay, cell_group = NULL, features_number_per_cell_group = 300){

  cell_group = enquo(cell_group)

  # Nest
  counts =
    counts |>
    nest(data = -!!cell_group)

  # If I have enough information per cell type
  if(counts |> filter(map_int(data, ncol) > 100) |> nrow() |> gt(1)){

    # Per cell type
    counts |>

      # Filter more than 10 cells
      filter(map_int(data, ncol) > 100) |>

      # Get feature within each cluster/cell-type
      mutate(feature = map(
        data,
        ~ .x |>
          FindVariableFeatures(nfeatures = features_number_per_cell_group, assay=assay) |>
          VariableFeatures(assay=assay)
      )) |>
      select(-data) |>
      unnest(feature) |>

      # Rename
      rename(group := !!cell_group)

  }
  else {

    warning(glue("jascap says: you have only one distinct `{quo_name(cell_group)}`, the per-cell-group variable gene detection will be skipped as it would olverlap with the global detection."))

    tibble()
  }
}

#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat VariableFeatures
seurat_to_variable_features_overall = function(counts, assay, features_number = 300){

  counts |>
    FindVariableFeatures(nfeatures = features_number, assay=assay) |>
    VariableFeatures(assay=assay) |>
    as_tibble() |>
    rename("feature" = "value") |>
    mutate(group= "variable_overall")

}

#' @importFrom rlang enquo
#' @importFrom rlang is_symbolic
#' @import tidyseurat
#' @importFrom Seurat NormalizeData
#'
#' @export
#'
#'
#'
seurat_to_variable_features = function(
    counts,
    assay,
    .sample,
    cell_group = NULL,
    features_number_independent_of_cell_groups = 300,
    features_number_per_cell_group = 300
){

  .sample = enquo(.sample)
  cell_group = enquo(cell_group)

  # If more than one sample balance the size
  if(counts |> distinct(!!.sample) |> nrow() |> gt(1))
    counts =
      counts |>

      # Sample up to a plateau to avoid extreme cell_type bias
      nest(data = -!!.sample) %>%
      mutate(n = map_int(data, ~ ncol(.x))) %>%
      mutate(upper_quantile = quantile(n, 0.75) %>% as.integer()) %>%
      mutate(data = map2(
        data, upper_quantile,
        ~ sample_n(.x, min(ncol(.x), .y), replace = FALSE)
      )) %>%
      filter(n>1) %>%
      unnest(data)

  # Normalise before - https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#normalizing-the-data-1
  counts  = counts |> NormalizeData(assay=assay)

  variable_df_overall = seurat_to_variable_features_overall(
    counts,
    assay,
    features_number = features_number_independent_of_cell_groups
  )

  # If cell_type_column_for_subsetting == null calculate clusters\
  if(!is_symbolic(cell_group)){

    counts =
      counts |>
      FindVariableFeatures(nfeatures = number_features_overall, assay=assay)  |>
      NormalizeData(assay=assay) |>
      ScaleData(assay=assay) |>
      RunPCA(assay=assay) |>
      FindNeighbors(dims = 1:20) |>
      FindClusters(resolution = 0.5)

    cell_group = as.symbol("seurat_clusters")
  }

  variable_df_by_cell_type = seurat_to_variable_features_by_cell_type(
    counts,
    assay,
    cell_group = !!cell_group,
    features_number_per_cell_group = features_number_per_cell_group
  )

  variable_df_overall  |>
    bind_rows(variable_df_by_cell_type)

}




# Cell Chat replacement
computeCommunProbPathway = function (object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05)
{
  if (is.null(net)) {
    net <- object@net
  }
  if (is.null(pairLR.use)) {
    pairLR.use <- object@LR$LRsig
  }
  prob <- net$prob
  prob[net$pval > thresh] <- 0
  pathways <- unique(pairLR.use$pathway_name)
  group <- factor(pairLR.use$pathway_name, levels = pathways)
  prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum),
                         c(2, 3, 1))
  pathways.sig <- pathways[apply(prob.pathways, 3, sum) !=
                             0]
  prob.pathways.sig <- prob.pathways[, , pathways.sig, drop=FALSE]
  idx <- sort(apply(prob.pathways.sig, 3, sum), decreasing = TRUE,
              index.return = TRUE)$ix
  pathways.sig <- pathways.sig[idx]
  prob.pathways.sig <- prob.pathways.sig[, , idx]
  if (is.null(object)) {
    netP = list(pathways = pathways.sig, prob = prob.pathways.sig)
    return(netP)
  }
  else {
    object@netP$pathways <- pathways.sig
    object@netP$prob <- prob.pathways.sig
    return(object)
  }
}

cellchat_process_sample_signal = function (object, signaling = NULL, pattern = c("outgoing",
                                                                                 "incoming", "all"), slot.name = "netP", color.use = NULL,
                                           color.heatmap = "BuGn", title = NULL, width = 10, height = 8,
                                           font.size = 8, font.size.title = 10, cluster.rows = FALSE,
                                           cluster.cols = FALSE)
{
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  }
  else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  }
  else if (pattern == "all") {
    mat <- t(outgoing + incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  }
  else {
    title <- paste0(paste0(legend.name, " signaling patterns"),
                    " - ", title)
  }
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
  mat[mat == 0] <- NA

  mat %>%
    as_tibble(rownames = "gene") %>%
    gather(cell_type, value, -gene) %>%
    mutate(value = if_else(value %in% c(NaN, NA), 0, value))
}

cellchat_circle_plot = function(pathway, x, y, DB, joint){

  cellchat_diff_for_circle(pathway, x, y) %>%
    draw_cellchat_circle_plot(
      vertex.weight = as.numeric((table(x@idents) + table(y@idents))/2),
      title.name = paste(pathway, DB, "\n", select_genes_for_circle_plot(joint, pathway)),
      edge.width.max = 4
    )
}

cellchat_matrix_for_circle = function (object, signaling, signaling.name = NULL, color.use = NULL,
                                       vertex.receiver = NULL, sources.use = NULL, targets.use = NULL,
                                       top = 1, remove.isolate = FALSE, vertex.weight = NULL, vertex.weight.max = NULL,
                                       vertex.size.max = 15, weight.scale = TRUE, edge.weight.max = NULL,
                                       edge.width.max = 8, layout = c("hierarchy", "circle", "chord"),
                                       thresh = 0.05, from = NULL, to = NULL, bidirection = NULL,
                                       vertex.size = NULL, pt.title = 12, title.space = 6, vertex.label.cex = 0.8,
                                       group = NULL, cell.order = NULL, small.gap = 1, big.gap = 10,
                                       scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20,
                                       legend.pos.y = 20, ...) {

  if(object@LR$LRsig %>% filter(pathway_name == signaling) %>% nrow %>% magrittr::equals(0)) return(NULL)

  layout <- match.arg(layout)
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents))
  }
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig,
                       key = "pathway_name", matching.exact = T, pair.only = T)
  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net
  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval
  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name],
                                         3, sum) != 0]
  }
  else {
    pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) !=
                                     0]
  }
  if (length(pairLR.name.use) == 0) {
    return(NULL)
    #stop(paste0("There is no significant communication of ", 					signaling.name))
  }
  else {
    pairLR <- pairLR[pairLR.name.use, ]
  }
  nRow <- length(pairLR.name.use)
  prob <- prob[, , pairLR.name.use]
  pval <- pval[, , pairLR.name.use]
  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify = "array")
    pval <- replicate(1, pval, simplify = "array")
  }
  prob.sum <- apply(prob, c(1, 2), sum)

  prob.sum
}

cellchat_diff_for_circle = function(pathway, x, y){
  zero_matrix =
    pathway %>%
    when(
      (.) %in% x@netP$pathways ~ cellchat_matrix_for_circle(x, layout = "circle", signaling = .),
      (.) %in% y@netP$pathways ~ cellchat_matrix_for_circle(y, layout = "circle", signaling = .)
    ) %>%
    `-` (.,.)

  m1 = pathway %>%
    when(
      (.) %in% x@netP$pathways ~ cellchat_matrix_for_circle(x, layout = "circle", signaling = .),
      ~ zero_matrix
    )

  m2 = pathway %>%
    when(
      (.) %in% y@netP$pathways ~ cellchat_matrix_for_circle(y, layout = "circle", signaling = .),
      ~ zero_matrix
    )

  m2 - m1
}

draw_cellchat_circle_plot = function (net, color.use = NULL, title.name = NULL, sources.use = NULL,
                                      targets.use = NULL, remove.isolate = FALSE, top = 1, top_absolute = NULL, weight.scale = T,
                                      vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = 15,
                                      vertex.label.cex = 0.8, vertex.label.color = "black", edge.weight.max = NULL,
                                      edge.width.max = 8, alpha.edge = 0.6, label.edge = FALSE,
                                      edge.label.color = "black", edge.label.cex = 0.8, edge.curved = 0.2,
                                      shape = "circle", layout = in_circle(), margin = 0.2, vertex.size = NULL,
                                      arrow.width = 1, arrow.size = 0.2)
{
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  options(warn = -1)

  if(!is.null(top_absolute)) {
    thresh = top_absolute
    net[abs(net) < thresh] <- 0
  }

  thresh <- stats::quantile(as.numeric(net) %>% abs %>% .[.>0], probs = 1 - top)

  net[abs(net) < thresh] <- 0

  if(sum(net)==0) return(NULL)

  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]],
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    if(length(idx)>0){
      net <- net[-idx, ,drop=FALSE]
      net <- net[, -idx, drop=FALSE]
    }
  }
  g <- graph_from_adjacency_matrix(net, mode = "directed",
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max +
    5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, -atan(coords_scale[igraph::V(g),
                                                                             2]/coords_scale[igraph::V(g), 1]), pi - atan(coords_scale[igraph::V(g),
                                                                                                                                       2]/coords_scale[igraph::V(g), 1]))
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(abs(igraph::E(g)$weight))
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + abs(igraph::E(g)$weight)/edge.weight.max *
      edge.width.max
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * abs(igraph::E(g)$weight)
  }
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex

  igraph::E(g)$color =
    circlize::colorRamp2(seq(max(abs(igraph::E(g)$weight)), -max(abs(igraph::E(g)$weight)), length.out =11), RColorBrewer::brewer.pal(11, "RdBu"))(igraph::E(g)$weight) %>%
    grDevices::adjustcolor(alpha.edge)


  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[,
                                                                1])] <- loop.angle[edge.start[which(edge.start[,
                                                                                                               2] == edge.start[, 1]), 1]]
  }
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)),
                               direction = -1, start = 0)
  label.dist <- vertex.weight/max(vertex.weight) + 2
  plot(g, edge.curved = edge.curved, vertex.shape = shape,
       layout = coords_scale, margin = margin, vertex.label.dist = label.dist,
       vertex.label.degree = label.locs, vertex.label.family = "Helvetica",
       edge.label.family = "Helvetica")
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 0.8)
  }

  grab_grob() |> as_grob() |> wrap_elements()

}

select_genes_for_circle_plot = function(x, pathway){
  paste(
    c(
      x@data.signaling[rownames(x@data.signaling) %in% (CellChatDB.human$interaction %>% filter(pathway_name == pathway) %>% distinct(ligand) %>% pull(1)),, drop=F] %>% rowSums() %>% .[(.)>100] %>% names(),
      x@data.signaling[rownames(x@data.signaling) %in% (CellChatDB.human$interaction %>% filter(pathway_name == pathway) %>% distinct(receptor) %>% pull(1)),, drop=F] %>% rowSums() %>% .[(.)>100] %>% names()
    ) %>% unique(),
    collapse = ","
  )

}

get_table_for_cell_vs_axis_bubble_plot = function (object, sources.use = NULL, targets.use = NULL, signaling = NULL,
                                                   pairLR.use = NULL, color.heatmap = c("Spectral", "viridis"),
                                                   n.colors = 10, direction = -1, thresh = 0.05, comparison = NULL,
                                                   group = NULL, remove.isolate = FALSE, max.dataset = NULL,
                                                   min.dataset = NULL, min.quantile = 0, max.quantile = 1,
                                                   line.on = TRUE, line.size = 0.2, color.text.use = TRUE,
                                                   color.text = NULL, title.name = NULL, font.size = 10, font.size.title = 10,
                                                   show.legend = TRUE, grid.on = TRUE, color.grid = "grey90",
                                                   angle.x = 90, vjust.x = NULL, hjust.x = NULL, return.data = FALSE)
{

  # cells.level <- levels(object@idents)
  # source.use.numerical = which(cells.level == source.use)
  #
  color.heatmap <- match.arg(color.heatmap)
  if (is.list(object@net[[1]])) {
    message("Comparing communications on a merged object \n")
  }
  else {
    message("Comparing communications on a single object \n")
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle = c(0, 45, 90)
    hjust = c(0, 1, 1)
    vjust = c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }
  if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
    })
  }
  else {
    color.use <- color.heatmap
  }
  if (direction == -1) {
    color.use <- rev(color.use)
  }
  if (is.null(comparison)) {
    cells.level <- levels(object@idents)
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }

    # TRY CATCH
    df.net <- 	tryCatch(
      expr = {
        subsetCommunication(object, slot.name = "net",
                            sources.use = sources.use, targets.use = targets.use,
                            signaling = signaling, pairLR.use = pairLR.use,
                            thresh = thresh)
      },
      error = function(e){
        return(NULL)
      }
    )

    if(is.null(df.net)) return(NULL)

    df.net$source.target <- paste(df.net$source, df.net$target,
                                  sep = " -> ")
    source.target <- paste(rep(sources.use, each = length(targets.use)),
                           targets.use, sep = " -> ")
    source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
    if (length(source.target.isolate) > 0) {
      df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate),
                                             ncol = ncol(df.net)))
      colnames(df.net.isolate) <- colnames(df.net)
      df.net.isolate$source.target <- source.target.isolate
      df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
      df.net.isolate$pval <- 1
      a <- stringr::str_split(df.net.isolate$source.target,
                              " -> ", simplify = T)
      df.net.isolate$source <- as.character(a[, 1])
      df.net.isolate$target <- as.character(a[, 2])
      df.net <- rbind(df.net, df.net.isolate)
    }
    df.net$pval[df.net$pval > 0.05] = 1
    df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
    df.net$pval[df.net$pval <= 0.01] = 3
    df.net$prob[df.net$prob == 0] <- NA
    df.net$prob.original <- df.net$prob
    df.net$prob <- -1/log(df.net$prob)
    idx1 <- which(is.infinite(df.net$prob) | df.net$prob <
                    0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.net$prob, na.rm = T) *
                             1.1, max(df.net$prob, na.rm = T) * 1.5, length.out = length(idx1))
      position <- sort(prob.original[idx1], index.return = TRUE)$ix
      df.net$prob[idx1] <- values.assign[match(1:length(idx1),
                                               position)]
    }
    df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in%
                                                                  unique(df.net$source)])
    df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in%
                                                                  unique(df.net$target)])
    group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))),
                         levels(df.net$target), sep = " -> ")
    df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
    df.net <- with(df.net, df.net[order(interaction_name_2),
    ])
    df.net$interaction_name_2 <- factor(df.net$interaction_name_2,
                                        levels = unique(df.net$interaction_name_2))
    cells.order <- group.names
    df.net$source.target <- factor(df.net$source.target,
                                   levels = cells.order)
    df <- df.net
  }
  else {
    dataset.name <- names(object@net)
    df.net.all <- subsetCommunication(object, slot.name = "net",
                                      sources.use = sources.use, targets.use = targets.use,
                                      signaling = signaling, pairLR.use = pairLR.use,
                                      thresh = thresh)
    df.all <- data.frame()
    for (ii in 1:length(comparison)) {
      cells.level <- levels(object@idents[[comparison[ii]]])
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- df.net.all[[comparison[ii]]]
      df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
      df.net$source.target <- paste(df.net$source, df.net$target,
                                    sep = " -> ")
      source.target <- paste(rep(sources.use, each = length(targets.use)),
                             targets.use, sep = " -> ")
      source.target.isolate <- setdiff(source.target,
                                       unique(df.net$source.target))
      if (length(source.target.isolate) > 0) {
        df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate),
                                               ncol = ncol(df.net)))
        colnames(df.net.isolate) <- colnames(df.net)
        df.net.isolate$source.target <- source.target.isolate
        df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
        df.net.isolate$pval <- 1
        a <- stringr::str_split(df.net.isolate$source.target,
                                " -> ", simplify = T)
        df.net.isolate$source <- as.character(a[, 1])
        df.net.isolate$target <- as.character(a[, 2])
        df.net <- rbind(df.net, df.net.isolate)
      }
      df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in%
                                                                    unique(df.net$source)])
      df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in%
                                                                    unique(df.net$target)])
      group.names <- paste(rep(levels(df.net$source),
                               each = length(levels(df.net$target))), levels(df.net$target),
                           sep = " -> ")
      group.names0 <- group.names
      group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]],
                            ")")
      if (nrow(df.net) > 0) {
        df.net$pval[df.net$pval > 0.05] = 1
        df.net$pval[df.net$pval > 0.01 & df.net$pval <=
                      0.05] = 2
        df.net$pval[df.net$pval <= 0.01] = 3
        df.net$prob[df.net$prob == 0] <- NA
        df.net$prob.original <- df.net$prob
        df.net$prob <- -1/log(df.net$prob)
      }
      else {
        df.net <- as.data.frame(matrix(NA, nrow = length(group.names),
                                       ncol = 5))
        colnames(df.net) <- c("interaction_name_2",
                              "source.target", "prob", "pval", "prob.original")
        df.net$source.target <- group.names0
      }
      df.net$group.names <- as.character(df.net$source.target)
      df.net$source.target <- paste0(df.net$source.target,
                                     " (", dataset.name[comparison[ii]], ")")
      df.net$dataset <- dataset.name[comparison[ii]]
      df.all <- rbind(df.all, df.net)
    }
    if (nrow(df.all) == 0) {
      return(NULL)
      #stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }
    idx1 <- which(is.infinite(df.all$prob) | df.all$prob <
                    0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.all$prob, na.rm = T) *
                             1.1, max(df.all$prob, na.rm = T) * 1.5, length.out = length(idx1))
      position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
      df.all$prob[idx1] <- values.assign[match(1:length(idx1),
                                               position)]
    }
    df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]
    df <- df.all
    df <- with(df, df[order(interaction_name_2), ])
    df$interaction_name_2 <- factor(df$interaction_name_2,
                                    levels = unique(df$interaction_name_2))
    cells.order <- c()
    dataset.name.order <- c()
    for (i in 1:length(group.names0)) {
      for (j in 1:length(comparison)) {
        cells.order <- c(cells.order, paste0(group.names0[i],
                                             " (", dataset.name[comparison[j]], ")"))
        dataset.name.order <- c(dataset.name.order,
                                dataset.name[comparison[j]])
      }
    }
    df$source.target <- factor(df$source.target, levels = cells.order)
  }
  min.cutoff <- quantile(df$prob, min.quantile, na.rm = T)
  max.cutoff <- quantile(df$prob, max.quantile, na.rm = T)
  df$prob[df$prob < min.cutoff] <- min.cutoff
  df$prob[df$prob > max.cutoff] <- max.cutoff
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (!is.null(max.dataset)) {
    signaling <- as.character(unique(df$interaction_name_2))
    for (i in signaling) {
      df.i <- df[df$interaction_name_2 == i, , drop = FALSE]
      cell <- as.character(unique(df.i$group.names))
      for (j in cell) {
        df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
        values <- df.i.j$prob
        idx.max <- which(values == max(values, na.rm = T))
        idx.min <- which(values == min(values, na.rm = T))
        dataset.na <- c(df.i.j$dataset[is.na(values)],
                        setdiff(dataset.name[comparison], df.i.j$dataset))
        if (length(idx.max) > 0) {
          if (!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset])) {
            df.i.j$prob <- NA
          }
          else if ((idx.max != idx.min) & !is.null(min.dataset)) {
            if (!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset])) {
              df.i.j$prob <- NA
            }
            else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in%
                                                    dataset.na)) > 0) {
              df.i.j$prob <- NA
            }
          }
        }
        df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
      }
      df[df$interaction_name_2 == i, "prob"] <- df.i$prob
    }
  }
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (nrow(df) == 0) {
    return(NULL)
    #stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
  }
  df$interaction_name_2 <- factor(df$interaction_name_2, levels = unique(df$interaction_name_2))
  df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target),
                                                                    unique(df$source.target)))
  df
}

grab_grob <- function(){
  grid.echo()
  grid.grab()
}
