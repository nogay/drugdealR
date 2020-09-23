#' Check if distance matrix of a disease changes after removing a given drug
#'
#' @param disease
#' @param drug
#' @param gr
#' @param subnet
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' library(drugdealR)
#' measure_dists('Myeloid Leukemia, Chronic', 'Imatinib', gr, subnet)
#' }

measure_dists <- function(disease, drug, gr, subnet){
  gr_noDrug <- subnet[-which(subnet$Protein_A == drug),]
  gr_noDrug <- igraph::graph_from_data_frame(gr_noDrug)
  gr_noDrug <- igraph::as.undirected(gr_noDrug)
  gr_noDrug <- igraph::simplify(gr_noDrug, remove.multiple = TRUE, remove.loops = TRUE)
  disNetDists <- sapply(disease, function(x) igraph::distances(gr,v=x,to=disease))
  disNetDists <- as.numeric(disNetDists)
  disNetDists_noDrug <- sapply(disease, function(x) igraph::distances(gr_noDrug,v=x,to=disease))
  rownames(disNetDists_noDrug) <- disease
  disNetDists_noDrug <- as.numeric(disNetDists_noDrug)
  return(all.equal(disNetDists,disNetDists_noDrug))
}
