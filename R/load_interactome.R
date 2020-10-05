#' Loads the protein-protein interaction network
#'
#' @return graph
#' @export
#'
#' @examples
#' \dontrun{
#' library(drugdealR)
#' load_interactome()
#' }

load_interactome <- function(){
  utils::data("interactome_subset")

  INTERACTOME <- igraph::graph_from_data_frame(interactome_subset)
  INTERACTOME <- igraph::as.undirected(INTERACTOME) #convert to undirected
  INTERACTOME <- igraph::simplify(INTERACTOME, remove.multiple = TRUE, remove.loops = TRUE)

  return(INTERACTOME)
}
