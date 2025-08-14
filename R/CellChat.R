# More robust implementation that does not fail if no results

#' @importFrom CellChat searchPair
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

  # Fix GCHECK 
  pathway_name = NULL 
  
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

#' @importFrom tidyr gather
#' @importFrom dplyr if_else
cellchat_process_sample_signal = function (object, signaling = NULL, pattern = c("outgoing",
                                                                                 "incoming", "all"), slot.name = "netP", color.use = NULL,
                                           color.heatmap = "BuGn", title = NULL, width = 10, height = 8,
                                           font.size = 8, font.size.title = 10, cluster.rows = FALSE,
                                           cluster.cols = FALSE)
{
  # Fix GCHECK 
  cell_type = NULL 
  value = NULL 
  gene = NULL 
  
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

computeCommunProbPathway= function (object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05)
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

  # STEFANO FIX
  if(length(levels(group))==1){
    xx = apply(prob, c(1, 2), by, group, sum)
    prob.pathways = xx |> array(dim = c(nrow(xx), ncol(xx), 1), dimnames = list(rownames(xx), colnames(xx), levels(group)))
  }
  else
    prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum),
                           c(2, 3, 1))



  pathways.sig <- pathways[apply(prob.pathways, 3, sum) !=     0]
  prob.pathways.sig <- prob.pathways[, , pathways.sig, drop=FALSE]
  idx <- sort(apply(prob.pathways.sig, 3, sum), decreasing = TRUE,
              index.return = TRUE)$ix
  pathways.sig <- pathways.sig[idx]
  prob.pathways.sig <- prob.pathways.sig[, , idx, drop=FALSE]
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

#' @importFrom CellChat triMean computeExpr_LR computeExpr_coreceptor computeRegionDistance computeExpr_agonist
computeCommunProb = function (object, type = c("triMean", "truncatedMean", "thresholdedMean",
                                               "median"), trim = 0.1, LR.use = NULL, raw.use = TRUE, population.size = FALSE,
                              distance.use = TRUE, interaction.length = 200, scale.distance = 0.01,
                              k.min = 10, nboot = 100, seed.use = 1L, Kh = 0.5, n = 1)
{
  
  
  type <- match.arg(type)
  cat(type, "is used for calculating the average gene expression per cell group.",
      "\n")
  FunMean <- switch(type, triMean = triMean, truncatedMean = function(x) mean(x,
                                                                              trim = trim, na.rm = TRUE), median = function(x) median(x,
                                                                                                                                      na.rm = TRUE))
  if (raw.use) {
    data <- as.matrix(object@data.signaling)
  }
  else {
    data <- object@data.project
  }
  if (is.null(LR.use)) {
    pairLR.use <- object@LR$LRsig
  }
  else {
    pairLR.use <- LR.use
  }
  complex_input <- object@DB$complex
  cofactor_input <- object@DB$cofactor
  my.sapply <-  sapply
  ptm = Sys.time()
  pairLRsig <- pairLR.use
  group <- object@idents
  geneL <- as.character(pairLRsig$ligand)
  geneR <- as.character(pairLRsig$receptor)
  nLR <- nrow(pairLRsig)
  numCluster <- nlevels(group)
  if (numCluster != length(unique(group))) {
    stop("Please check `unique(object@idents)` and ensure that the factor levels are correct!\n         You may need to drop unused levels using 'droplevels' function. e.g.,\n         `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`")
  }
  data.use <- data/max(data)
  nC <- ncol(data.use)
  data.use.avg <- aggregate(t(data.use), list(group), FUN = FunMean)
  data.use.avg <- t(data.use.avg[, -1])
  colnames(data.use.avg) <- levels(group)
  dataLavg <- computeExpr_LR(geneL, data.use.avg, complex_input)
  dataRavg <- computeExpr_LR(geneR, data.use.avg, complex_input)

  # STEFANO ADDED THIS
  dataRavg[dataRavg=="NaN"]=0
  dataLavg[dataLavg=="NaN"]=0

  dataRavg.co.A.receptor <- computeExpr_coreceptor(cofactor_input,
                                                   data.use.avg, pairLRsig, type = "A")
  dataRavg.co.I.receptor <- computeExpr_coreceptor(cofactor_input,
                                                   data.use.avg, pairLRsig, type = "I")
  dataRavg <- dataRavg * dataRavg.co.A.receptor/dataRavg.co.I.receptor
  dataLavg2 <- t(replicate(nrow(dataLavg), as.numeric(table(group))/nC))
  dataRavg2 <- dataLavg2
  index.agonist <- which(!is.na(pairLRsig$agonist) & pairLRsig$agonist !=
                           "")
  index.antagonist <- which(!is.na(pairLRsig$antagonist) &
                              pairLRsig$antagonist != "")
  if (object@options$datatype != "RNA") {
    data.spatial <- object@images$coordinates
    spot.size.fullres <- object@images$scale.factors$spot
    spot.size <- object@images$scale.factors$spot.diameter
    d.spatial <- computeRegionDistance(coordinates = data.spatial,
                                       group = group, trim = trim, interaction.length = interaction.length,
                                       spot.size = spot.size, spot.size.fullres = spot.size.fullres,
                                       k.min = k.min)
    if (distance.use) {
      print(paste0(">>> Run CellChat on spatial imaging data using distances as constraints <<< [",
                   Sys.time(), "]"))
      d.spatial <- d.spatial * scale.distance
      diag(d.spatial) <- NaN
      cat("The suggested minimum value of scaled distances is in [1,2], and the calculated value here is ",
          min(d.spatial, na.rm = TRUE), "\n")
      if (min(d.spatial, na.rm = TRUE) < 1) {
        stop("Please increase the value of `scale.distance` and check the suggested values in the parameter description (e.g., 1, 0.1, 0.01, 0.001, 0.11, 0.011)")
      }
      P.spatial <- 1/d.spatial
      P.spatial[is.na(d.spatial)] <- 0
      diag(P.spatial) <- max(P.spatial)
      d.spatial <- d.spatial/scale.distance
    }
    else {
      print(paste0(">>> Run CellChat on spatial imaging data without distances as constraints <<< [",
                   Sys.time(), "]"))
      P.spatial <- matrix(1, nrow = numCluster, ncol = numCluster)
      P.spatial[is.na(d.spatial)] <- 0
    }
  }
  else {
    print(paste0(">>> Run CellChat on sc/snRNA-seq data <<< [",
                 Sys.time(), "]"))
    d.spatial <- matrix(NaN, nrow = numCluster, ncol = numCluster)
    P.spatial <- matrix(1, nrow = numCluster, ncol = numCluster)
    distance.use = NULL
    interaction.length = NULL
    spot.size = NULL
    spot.size.fullres = NULL
    k.min = NULL
  }
  Prob <- array(0, dim = c(numCluster, numCluster, nLR))
  Pval <- array(0, dim = c(numCluster, numCluster, nLR))
  set.seed(seed.use)
  permutation <- replicate(nboot, sample.int(nC, size = nC))
  data.use.avg.boot <- my.sapply(X = 1:nboot, FUN = function(nE) {
    groupboot <- group[permutation[, nE]]
    data.use.avgB <- aggregate(t(data.use), list(groupboot),
                               FUN = FunMean)
    data.use.avgB <- t(data.use.avgB[, -1])
    return(data.use.avgB)
  }, simplify = FALSE)
  pb <- txtProgressBar(min = 0, max = nLR, style = 3, file = stderr())
  for (i in 1:nLR) {
    dataLR <- Matrix::crossprod(matrix(dataLavg[i, ], nrow = 1),
                                matrix(dataRavg[i, ], nrow = 1))
    P1 <- dataLR^n/(Kh^n + dataLR^n)
    P1_Pspatial <- P1 * P.spatial
    if (sum(P1_Pspatial) == 0) {
      Pnull = P1_Pspatial
      Prob[, , i] <- Pnull
      p = 1
      Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster,
                            byrow = FALSE)
    }
    else {
      if (is.element(i, index.agonist)) {
        data.agonist <- computeExpr_agonist(data.use = data.use.avg,
                                            pairLRsig, cofactor_input, index.agonist = i,
                                            Kh = Kh, n = n)
        P2 <- Matrix::crossprod(matrix(data.agonist,
                                       nrow = 1))
      }
      else {
        P2 <- matrix(1, nrow = numCluster, ncol = numCluster)
      }
      if (is.element(i, index.antagonist)) {
        data.antagonist <- CellChat::computeExpr_antagonist(data.use = data.use.avg,
                                                  pairLRsig, cofactor_input, index.antagonist = i,
                                                  Kh = Kh, n = n)
        P3 <- Matrix::crossprod(matrix(data.antagonist,
                                       nrow = 1))
      }
      else {
        P3 <- matrix(1, nrow = numCluster, ncol = numCluster)
      }
      if (population.size) {
        P4 <- Matrix::crossprod(matrix(dataLavg2[i,
        ], nrow = 1), matrix(dataRavg2[i, ], nrow = 1))
      }
      else {
        P4 <- matrix(1, nrow = numCluster, ncol = numCluster)
      }
      Pnull = P1 * P2 * P3 * P4 * P.spatial
      Prob[, , i] <- Pnull
      Pnull <- as.vector(Pnull)
      Pboot <- sapply(X = 1:nboot, FUN = function(nE) {
        data.use.avgB <- data.use.avg.boot[[nE]]
        dataLavgB <- computeExpr_LR(geneL[i], data.use.avgB,
                                    complex_input)
        dataRavgB <- computeExpr_LR(geneR[i], data.use.avgB,
                                    complex_input)
        dataRavgB.co.A.receptor <- computeExpr_coreceptor(cofactor_input,
                                                          data.use.avgB, pairLRsig[i, , drop = FALSE],
                                                          type = "A")
        dataRavgB.co.I.receptor <- computeExpr_coreceptor(cofactor_input,
                                                          data.use.avgB, pairLRsig[i, , drop = FALSE],
                                                          type = "I")
        dataRavgB <- dataRavgB * dataRavgB.co.A.receptor/dataRavgB.co.I.receptor
        dataLRB = Matrix::crossprod(dataLavgB, dataRavgB)
        P1.boot <- dataLRB^n/(Kh^n + dataLRB^n)
        if (is.element(i, index.agonist)) {
          data.agonist <- computeExpr_agonist(data.use = data.use.avgB,
                                              pairLRsig, cofactor_input, index.agonist = i,
                                              Kh = Kh, n = n)
          P2.boot <- Matrix::crossprod(matrix(data.agonist,
                                              nrow = 1))
        }
        else {
          P2.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
        }
        if (is.element(i, index.antagonist)) {
          data.antagonist <- CellChat::computeExpr_antagonist(data.use = data.use.avgB,
                                                    pairLRsig, cofactor_input, index.antagonist = i,
                                                    Kh = Kh, n = n)
          P3.boot <- Matrix::crossprod(matrix(data.antagonist,
                                              nrow = 1))
        }
        else {
          P3.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
        }
        if (population.size) {
          groupboot <- group[permutation[, nE]]
          dataLavg2B <- as.numeric(table(groupboot))/nC
          dataLavg2B <- matrix(dataLavg2B, nrow = 1)
          dataRavg2B <- dataLavg2B
          P4.boot = Matrix::crossprod(dataLavg2B, dataRavg2B)
        }
        else {
          P4.boot = matrix(1, nrow = numCluster, ncol = numCluster)
        }
        Pboot = P1.boot * P2.boot * P3.boot * P4.boot *
          P.spatial
        return(as.vector(Pboot))
      })
      Pboot <- matrix(unlist(Pboot), nrow = length(Pnull),
                      ncol = nboot, byrow = FALSE)
      nReject <- rowSums(Pboot - Pnull > 0)
      p = nReject/nboot
      Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster,
                            byrow = FALSE)
    }
    setTxtProgressBar(pb = pb, value = i)
  }
  close(con = pb)
  Pval[Prob == 0] <- 1
  dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
  dimnames(Pval) <- dimnames(Prob)
  net <- list(prob = Prob, pval = Pval)
  execution.time = Sys.time() - ptm
  object@options$run.time <- as.numeric(execution.time, units = "secs")
  object@options$parameter <- list(type.mean = type, trim = trim,
                                   raw.use = raw.use, population.size = population.size,
                                   nboot = nboot, seed.use = seed.use, Kh = Kh, n = n,
                                   distance.use = distance.use, interaction.length = interaction.length,
                                   spot.size = spot.size, spot.size.fullres = spot.size.fullres,
                                   k.min = k.min)
  if (object@options$datatype != "RNA") {
    object@images$distance <- d.spatial
  }
  object@net <- net
  print(paste0(">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [",
               Sys.time(), "]"))
  return(object)
}

#' @importFrom future nbrOfWorkers
#' @importFrom methods slot
#' @importFrom pbapply pbsapply
#' @import future.apply
netAnalysis_computeCentrality = function (object = NULL, slot.name = "netP", net = NULL, net.name = NULL,
                                          thresh = 0.05)
{
  if (is.null(net)) {
    prob <- methods::slot(object, slot.name)$prob
    pval <- methods::slot(object, slot.name)$pval
    pval[prob == 0] <- 1
    prob[pval >= thresh] <- 0
    net = prob
  }
  if (is.null(net.name)) {
    net.name <- dimnames(net)[[3]]
  }
  if (length(dim(net)) == 3) {
    nrun <- dim(net)[3]
    my.sapply <- ifelse(test = future::nbrOfWorkers() ==
                          1, yes = pbapply::pbsapply, no = future.apply::future_sapply)
    centr.all = my.sapply(X = 1:nrun, FUN = function(x) {
      net0 <- net[, , x]

      # ADDED BY STEFANO
      net0[net0<0] = 0

      return(computeCentralityLocal(net0))
    }, simplify = FALSE)
  }
  else {
    centr.all <- as.list(computeCentralityLocal(net))
  }
  names(centr.all) <- net.name
  if (is.null(object)) {
    return(centr.all)
  }
  else {
    slot(object, slot.name)$centr <- centr.all
    return(object)
  }
}

cellchat_circle_plot = function(pathway, x, y, DB, joint){

  cellchat_diff_for_circle(pathway, x, y) %>%
    draw_cellchat_circle_plot(
      vertex.weight = as.numeric((table(x@idents) + table(y@idents))/2),
      title.name = paste(pathway, DB, "\n", select_genes_for_circle_plot(joint, pathway)),
      edge.width.max = 4
    )
}

#' @importFrom purrr when
cellchat_diff_for_circle = function(pathway, x, y){
  
  # Fix GCHECK 
  . = NULL
  
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

#' @importFrom igraph graph_from_adjacency_matrix layout_
#' @importFrom CellChat scPalette
#' @importFrom reshape2 melt
#' @importFrom patchwork wrap_elements
#' @importFrom cowplot as_grob
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales rescale
#' @importFrom igraph in_circle 
draw_cellchat_circle_plot = function (net, color.use = NULL, title.name = NULL, sources.use = NULL,
                                      targets.use = NULL, remove.isolate = FALSE, top = 1, top_absolute = NULL, weight.scale = T,
                                      vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = 15,
                                      vertex.label.cex = 0.8, vertex.label.color = "black", edge.weight.max = NULL,
                                      edge.width.max = 8, alpha.edge = 0.6, label.edge = FALSE,
                                      edge.label.color = "black", edge.label.cex = 0.8, edge.curved = 0.2,
                                      shape = "circle", layout = in_circle(), margin = 0.2, vertex.size = NULL,
                                      arrow.width = 1, arrow.size = 0.2)
{
  # Pass GCHECKS 
  target = NULL 
  
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
  g <- igraph::graph_from_adjacency_matrix(net, mode = "directed",
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- igraph::layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = CellChat::scPalette(length(igraph::V(g)))
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

  grab_grob() |> cowplot::as_grob() |> patchwork::wrap_elements()

}

#' @importFrom dplyr distinct
#' 
select_genes_for_circle_plot = function(x, pathway){
  # Fix GChecks 
  CellChatDB.human <- NULL
  pathway_name <- NULL
  ligand <- NULL
  . <- NULL
  receptor <- NULL
  

  paste(
    c(
      x@data.signaling[rownames(x@data.signaling) %in% (CellChatDB.human$interaction %>% filter(pathway_name == pathway) %>% distinct(ligand) %>% pull(1)),, drop=F] %>% rowSums() %>% .[(.)>100] %>% names(),
      x@data.signaling[rownames(x@data.signaling) %in% (CellChatDB.human$interaction %>% filter(pathway_name == pathway) %>% distinct(receptor) %>% pull(1)),, drop=F] %>% rowSums() %>% .[(.)>100] %>% names()
    ) %>% unique(),
    collapse = ","
  )

}

#' @importFrom CellChat subsetCommunication
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales viridis_pal
#' 
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
  
  # Fix GChecks 
  prob.original = NULL 

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
        CellChat::subsetCommunication(object, slot.name = "net",
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
      df.net.isolate <- BiocGenerics::as.data.frame(matrix(NA, nrow = length(source.target.isolate),
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
    df.net.all <- CellChat::subsetCommunication(object, slot.name = "net",
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
        df.net.isolate <- BiocGenerics::as.data.frame(matrix(NA, nrow = length(source.target.isolate),
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
        df.net <- BiocGenerics::as.data.frame(matrix(NA, nrow = length(group.names),
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

#' @importFrom gridGraphics grid.echo
#' @importFrom grid grid.grab
#' 
grab_grob <- function(){
  grid.echo()
  grid.grab()
}


#' Ligand-Receptor Count from Seurat Data
#'
#' @description
#' Calculates ligand-receptor interactions for each cell type in a Seurat object using CellChat.
#'
#' @param counts Seurat object.
#' @param .cell_group Cell group variable.
#' @param assay Name of the assay to use.
#' @param sample_for_plotting Sample name for plotting.
#'
#' @return A list of communication results including interactions and signaling pathways.
#'
#' @importFrom CellChat createCellChat
#' @importFrom CellChat setIdent
#' @importFrom CellChat subsetDB
#' @importFrom CellChat subsetData
#' @importFrom CellChat identifyOverExpressedGenes
#' @importFrom CellChat identifyOverExpressedInteractions
#' @importFrom CellChat smoothData
#' @importFrom CellChat filterCommunication
#' @importFrom CellChat aggregateNet
#' @importFrom rlang quo_name
#' @importFrom rlang enquo
#' @importFrom tibble tibble
#' @importFrom purrr map2
#' @importFrom purrr map
#' @importFrom purrr map2_dbl
#' @importFrom dplyr distinct add_count
#' @export
seurat_to_ligand_receptor_count = function(counts, .cell_group, assay, sample_for_plotting = ""){
  
  #Fix GChecks 
  cell_type_harmonised <- NULL
  n_cells <- NULL
  DB <- NULL
  cell_vs_all_cells_per_pathway <- NULL
  gene <- NULL
  
  # Your code for seurat_to_ligand_receptor_count function here
  
  
  .cell_group = enquo(.cell_group)
  
  # If only one cell, return empty
  if((counts |> distinct(!!.cell_group) |> nrow()) < 2) return(tibble)
  
  counts_cellchat =
    counts |>
    
    # Filter
    filter(!is.na(!!.cell_group)) |>
    
    # Filter cell types with > 10 cells
    add_count(cell_type_harmonised, name = "n_cells") |>
    filter(n_cells>=10) |>
    
    
    # Convert from seurat MUST BE LOG-NORMALISED
    createCellChat(group.by = quo_name(.cell_group), assay = assay) |>
    setIdent( ident.use = quo_name(.cell_group))
  
  
  communication_results =
    tibble(DB = c("Secreted Signaling", "ECM-Receptor" , "Cell-Cell Contact" )) |>
    mutate(data = list(counts_cellchat)) |>
    mutate(data = map2(
      data, DB,
      ~ {
        print(.y)
        .x@DB <- subsetDB(CellChat::CellChatDB.human, search = .y)
        
        x = .x |>
          subsetData() |>
          identifyOverExpressedGenes() |>
          identifyOverExpressedInteractions() |>
          smoothData(CellChat::PPI.human)
        
        if(nrow(x@LR$LRsig)==0) return(NA)
        
        x |>
          computeCommunProb() |>
          filterCommunication() |>
          computeCommunProbPathway() |>
          aggregateNet()
        
      }
    )) |>
    
    # Record sample
    mutate(sample = sample_for_plotting) |>
    
    # Add histogram
    mutate(tot_interactions = map2_dbl(
      data, DB,
      ~  .x |> when(
        !is.na(.) ~ sum(.x@net$count),
        ~ 0
      )
    )) |>
    
    # Add histogram
    mutate(cell_cell_count = map2_dbl(
      data, DB,
      ~  {
        my_data = .x
        
        # Return empty if no results
        if(is.na(my_data)) return(tibble(cell_from = character(), cell_to = character(),  weight = numeric()))
        
        my_data@net$count |>
          as_tibble(rownames = "cell_from") |>
          pivot_longer(-cell_from, names_to = "cell_to", values_to = "count")
        
        
      }
    )) |>
    
    # values_df_for_heatmap
    # Scores for each cell types across all others. How communicative is each cell type
    mutate(cell_vs_all_cells_per_pathway = map2(
      data  , sample,
      ~ when(
        .x,
        !is.na(.x) && length(.x@netP$pathways) > 0 ~
          netAnalysis_computeCentrality(., slot.name = "netP") |>
          cellchat_process_sample_signal(
            pattern = "all", signaling = .x@netP$pathways,
            title = .y, width = 5, height = 6, color.heatmap = "OrRd"
          ),
        ~ tibble(gene = character(), cell_type = character(), value = double())
      )
    ))
  
  genes = communication_results |> select(cell_vs_all_cells_per_pathway) |> unnest(cell_vs_all_cells_per_pathway) |> distinct(gene) |> pull(gene)
  
  # Hugh resolution
  communication_results |>
    mutate(cell_vs_cell_per_pathway = map(
      data,
      ~ {
        my_data = .x
        
        # Return empty if no results
        if(is.na(my_data)) return(tibble(gene = character(),  result = list()))
        
        tibble(gene = genes) |>
          mutate(result = map(gene, ~ {
            
            unparsed_result = cellchat_matrix_for_circle(my_data,  layout = "circle", signaling = .x)
            
            if(!is.null(unparsed_result))
              unparsed_result |>
              as_tibble(rownames = "cell_type_from") |>
              pivot_longer(-cell_type_from, names_to = "cell_type_to", values_to = "score")
            
            unparsed_result
            
          }))
      }
    )) |>
    
    mutate(cell_cell_weight = map(
      data,
      ~ {
        my_data = .x
        
        # Return empty if no results
        if(is.na(my_data)) return(tibble(cell_from = character(), cell_to = character(),  weight = numeric()))
        
        my_data@net$weight |>
          as_tibble(rownames = "cell_from") |>
          pivot_longer(-cell_from, names_to = "cell_to", values_to = "weight")
        
        
      }
    ))
  
  
}
