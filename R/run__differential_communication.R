
# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_path_ligand_receptor_count_output = args[2:(length(args)-5)]
input_metadata = args[[length(args)-4]]
output_path_plot_overall = args[[length(args)-3]]
output_path_plot_heatmap = args[[length(args)-2]]
output_path_plot_circle = args[[length(args)-1]]
output_path_values_communication = args[[length(args)]]

renv::load(project = code_directory)

library(tidyverse)
library(Seurat)
library(glue)
library(CellChat)
library(tidyseurat)
library(tidySingleCellExperiment)
library(tidyHeatmap)
library(purrr)
library(patchwork)
library(cowplot)
library(gridGraphics)

source(glue("{code_directory}/R/in_house_functions.R"))

# Create dir
output_path_plot_overall |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

factor_of_interest = as.symbol("factor_of_interest")

# input_path_ligand_receptor_count_output = dir("data/fast_pipeline_results/communication", full.names = TRUE)

# factor_of_interest = "severity"

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

integrated_cellchat_plot_df =
  input_path_ligand_receptor_count_output |>
  map_df(~ readRDS(.x)) |>

  # Add metadata
  left_join( readRDS(input_metadata), by = "sample") |>

  # Filter if not results otherwise it fails
  filter(map_lgl(data, ~ !is.na(.x))) |>

  # Nest
  nest(data_DB = -DB) |>

  mutate(joint = map(data_DB, ~ {
    object.list <- .x$data
    mergeCellChat(object.list, add.names = paste(.x$sample, pull(.x, severity), sep="___"))
  })) |>

  # Add histogram
  mutate(plot_histogram = map2(joint, DB, ~ {
   s =  sapply(.x@net, function(x) sum(x$count))

    when(is.na(s) ~ 0 , ~s) |>
    as.data.frame() |>
    as_tibble(rownames = "severity") |>
    setNames(c( "severity", "tot"))
  })) |>
  select(-joint) |>

  # Information flow
  # Add histogram
  mutate(plot_information_flow = map2(plot_histogram,DB, ~ {

    .x	|>
      separate(severity, c(  "sample", "severity"), "___") |>
      ggplot(aes(severity, tot, fill=severity)) +
      geom_boxplot() +
      geom_point() +
      ggpubr::stat_compare_means() +
      ggtitle(.y)
  }
  ))

# Plot overall counts
plot_overall_counts =
  integrated_cellchat_plot_df |>
  pull(plot_information_flow) |>
  wrap_plots(nrow = 1) +
  plot_layout(guides = 'collect' ) &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'), legend.position="bottom")



netAnalysis_computeCentrality =
  function (object = NULL, slot.name = "netP", net = NULL, net.name = NULL,  thresh = 0.05){
  if (is.null(net)) {
    prob <- methods::slot(object, slot.name)$prob
    pval <- methods::slot(object, slot.name)$pval
    pval[prob == 0] <- 1
    prob[pval >= thresh] <- 0
    net = prob
  }

    # If result is only one
    if(length(dimnames(net))<3)
      net = abind::abind(net, along = 3)

  if (is.null(net.name)) {
    net.name <- dimnames(net)[[3]]
  }
  if (length(dim(net)) == 3) {
    nrun <- dim(net)[3]
    my.sapply <- ifelse(test = future::nbrOfWorkers() ==
                          1, yes = pbapply::pbsapply, no = future.apply::future_sapply)
    centr.all = my.sapply(X = 1:nrun, FUN = function(x) {
      net0 <- net[, , x]
      return(CellChat:::computeCentralityLocal(net0))
    }, simplify = FALSE)
  }
  else {
    centr.all <- as.list(CellChat:::computeCentralityLocal(net))
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


values_df_for_heatmap =
  integrated_cellchat_plot_df |>
  select(DB, data_DB) |>
  unnest(data_DB) |>
  mutate(data = map2(
    data  , sample,
    ~ when(
      .x,
      length(.x@netP$pathways) > 0 ~
        netAnalysis_computeCentrality(., slot.name = "netP") |>
        cellchat_process_sample_signal(
          pattern = "all", signaling = .x@netP$pathways,
          title = .y, width = 5, height = 6, color.heatmap = "OrRd"
        ),
      ~ tibble(gene = character(), cell_type = character(), value = double())
      )
  )) |>
  unnest(data)  |>
  complete(nesting(DB, gene), nesting(severity, sample), cell_type, fill = list(value = 0))


signal_df_for_heatmap =
  values_df_for_heatmap |>
  nest(value_data = -c(DB  ,   gene ,  severity,  cell_type)) |>
  mutate(value = map_dbl(value_data, ~ mean(.x$value))) |>

  select(-value_data ) |>
  spread(severity, value) |>

  # THIS SHOULD CHANGE WITH THE RIGHT VALUES OF FACTOR OF INTEREST
  mutate(diff = severe - moderate) |>

  # Fix CD40
  mutate(gene = if_else(gene=="CD40" & DB=="Cell-Cell Contact" , "CD40_",  gene)) |>
  group_by(gene) |>
  mutate(max_diff = max(abs(diff))) |>
  filter(abs(max_diff)>0.34)

# Plot
plot_heatmap =
  signal_df_for_heatmap |>
  group_by(cell_type) |>
  mutate(sum_cell_diff = sum(diff)) |>
  group_by(gene) |>
  mutate(sum_gene_diff = sum(diff)) |>
  group_by(DB) |>

   mutate(gene = if_else(gene=="CD40" & DB=="Cell-Cell Contact" , "CD40_",  gene)) |>

  heatmap(
    gene, cell_type, diff,
    palette_value = circlize::colorRamp2(seq(0.5, -0.5, length.out =11), RColorBrewer::brewer.pal(11, "RdBu")),
    scale="none"
  ) |>
  add_bar(sum_cell_diff) |>
  add_bar(sum_gene_diff)

select_genes_for_circle_plot = function(x, pathway){
  paste(
    c(
      x@data.signaling[rownames(x@data.signaling) %in% (CellChatDB$interaction %>% filter(pathway_name == pathway) %>% distinct(ligand) %>% pull(1)),, drop=F] %>% rowSums() %>% .[(.)>100] %>% names(),
      x@data.signaling[rownames(x@data.signaling) %in% (CellChatDB$interaction %>% filter(pathway_name == pathway) %>% distinct(receptor) %>% pull(1)),, drop=F] %>% rowSums() %>% .[(.)>100] %>% names()
    ) %>% unique(),
    collapse = ","
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

CellChatDB <- CellChatDB.human

# plot_circle_communication_all
plot_circle_communication_all =
  integrated_cellchat_plot_df |>

  unnest(data_DB) |>
  left_join(signal_df_for_heatmap |> distinct(DB, gene) |> ungroup()) |>
  mutate(genes_in_pathway = map2(gene, data, ~ select_genes_for_circle_plot(.y, .x))) |>

  mutate(plot_cell_Cell_pathway = map2( data, gene , ~
                                          cellchat_matrix_for_circle(.x,  layout = "circle", signaling = .y)
  ))  |>
  filter(map_lgl(plot_cell_Cell_pathway, ~ !is.null(.x))) |>

  # Create the full column name because cellchat sucks
  mutate(all_cell_types = map(plot_cell_Cell_pathway, ~ colnames(.x))) |>
  mutate(all_cell_types = unlist(all_cell_types) |> unique() |> list()) |>
  mutate(missing_column = map2(
    plot_cell_Cell_pathway,
    all_cell_types,
    ~ .y |> setdiff(colnames(.x))
  )) |>
  mutate(plot_cell_Cell_pathway =
           map2(
             plot_cell_Cell_pathway,
             missing_column,
             ~ when(.y,
                    length(.) == 0 ~ .x,
                    ~ {
                      x = as.data.frame(.x)
                      x[.y] = 0
                      x[.y, ] = 0
                      x = as.matrix(x)

                    })
           )) |>
  mutate(plot_cell_Cell_pathway = map2(plot_cell_Cell_pathway, all_cell_types, ~ .x[.y, .y])) |>


  # Add line weights
  mutate(
    line_weights = map2(
      data,
      all_cell_types,
      ~ .y |>
        purrr::map_dfc(setNames, object = list(integer())) |>
        bind_rows(
          .x@idents |> table() |>  enframe() |> spread(name, value) |> mutate(across(everything(), ~ as.integer(.x)))
        ) |>
        mutate(across(everything(), ~ replace_na(.x, 0)))
    )
  ) |>

  # Summarise
  nest(data_path_type = -c(gene, severity)) |>

  mutate(
    pathway_mean = map(data_path_type, ~ purrr::reduce(.x$plot_cell_Cell_pathway, `+`)),
    line_weights_sum = map(
      data_path_type,
      ~ .x$line_weights |> purrr::reduce(bind_rows) |> colSums()
    ),
    genes_in_pathway = map(
      data_path_type,
      ~ .x$genes_in_pathway |> unlist() |> unique()
    )
  ) |>
  arrange(severity) |>

  # JUST SELECT COVID PATIENTS FOR NOW FOR TESTING
  filter(severity != "NA") |>

  # ADAPT TO OLD CODE
  mutate(severity = as.character(severity)) |>

  nest(data_type_path = -c(gene)) |>
  mutate(
    plot_diff = map(
      data_type_path,
      ~  when(
        pull(.x, severity),
        unlist(.) |>  identical("moderate") ~ -.x$pathway_mean[[1]],
        unlist(.) |>  identical("severe") ~ .x$pathway_mean[[1]],
        unlist(.) |> identical(c("moderate", "severe")) ~ .x$pathway_mean[[2]] - .x$pathway_mean[[1]],
        ~ stop("something went wrong")
      )
    ),
    line_weights_sum_sum = map(
      data_type_path,
      ~ .x$line_weights_sum |> purrr::reduce(`+`) |> magrittr::divide_by(2)
    ),
    genes_in_pathway = map(
      data_type_path,
      ~ .x$genes_in_pathway |> unlist() |> unique()
    )
  ) |>
  mutate(plot_diff_max = map_dbl(plot_diff, ~ max(abs(.x)))) |>
  mutate(plot_diff_quant = quantile(plot_diff_max, 0.25))






# Save
plot_overall_counts |> saveRDS(output_path_plot_overall)
plot_heatmap |> saveRDS(output_path_plot_heatmap)
plot_circle_communication_all |> saveRDS(output_path_plot_circle)
list(
  signal_df_for_heatmap = signal_df_for_heatmap,
  integrated_cellchat_plot_df = integrated_cellchat_plot_df
) |>
  saveRDS(output_path_values_communication)


# ggsave(
#   plot = purrr::reduce(p, `+`) + plot_layout( ncol = ceiling(sqrt(length(p))), nrow = ceiling(sqrt(length(p)))),
#   filename = output_path_plot_circle,
#   units = "mm",
#   width = 60*sqrt(length(p)),
#   height = 60*sqrt(length(p)),
#   device = "pdf"
# )


