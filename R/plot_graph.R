#' plot_graph
#' @description Simply function for visualizing the estimated network structure
#'
#' @param pcor_mat partial correlation matrix from GGMprojpred
#' @param pos_col  color for positive values
#' @param neg_col  color for negative values
#' @param scl      controls scale of edge sizes
#' @param ...      allows for passing agruments to plot.igraph
#'
#' @export
#'
#' @examples
plot_graph <- function(pcor_mat, pos_col = "green", neg_col = "red", scl = 10, ...){
  diag(pcor_mat) <- 0

  graph <-graph.adjacency(pcor_mat,weighted = TRUE, mode="undirected" )

  E(graph)$color <- ifelse(E(graph)$weight > 0, pos_col, neg_col)

  plot.igraph(graph, edge.width = abs(E(graph)$weight) * scl,
              edge.color = E(graph)$color,...)

}
