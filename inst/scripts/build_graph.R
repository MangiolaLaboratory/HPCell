library(igraph)

# load knowledge graph
adj_immune = read.csv("inst/extdata/immune_tree.csv", row.names = 1, check.names = FALSE) |>
  as.matrix()
immune_graph = graph_from_adjacency_matrix(adj_immune, mode = "directed", weighted = TRUE)

# Hierarchy - true (known) hierarchical relationships
# Consensus-only - temporary relationship used for ambiguous or often confused classes (e.g. macrophages and monocytes are often considered a single class but macrophages are not monocytes)
E(immune_graph)$Type = ifelse(E(immune_graph)$weight == 1, "Hierarchy", "Consensus-only")
E(immune_graph)$weight = 1
usethis::use_data(immune_graph)
