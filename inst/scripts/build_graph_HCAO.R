library(igraph)

# load knowledge graph
HCAO_edges = read.csv("inst/extdata/HCAO_immune_graph_edges.csv",stringsAsFactors = FALSE)
HCAO_vertices = read.csv("inst/extdata/HCAO_immune_graph_vertices.csv",stringsAsFactors = FALSE)
HCAO_graph = graph_from_data_frame(d = HCAO_edges, vertices = HCAO_vertices, directed = TRUE)

usethis::use_data(HCAO_graph)
