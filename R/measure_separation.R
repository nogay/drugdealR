#' Measure separation of two entities (disease or drugs)
#'
#' @param genesA list of genes associated with the given drug/disease.
#' @param genesB list of genes associated with the given drug/disease.
#' @param gr the subnetwork igraph graph object that contains all genes involved in
#' diseases (according to DisNor) and drug targets (according to Drugbank) of
#' interest.
#'
#' @return numeric
#' @export
#'
#' @examples
#' \dontrun{
#' library(drugdealR)
#' d <- load_interactome()
#'
#' To calculate separation between two drugs, insert drug names
#' measure_separation(d$drug1_genes, d$drug2_genes, gr)
#'
#' To calculate separation between two diseases, insert disease names
#' measure_separation(d$disease1_genes, d$disease2_genes, gr)
#' }

measure_separation <- function(genesA, genesB, gr){
  dAA <- sapply(genesA, function(x) igraph::distances(gr,v=x,to=genesA))
  dAA[dAA == Inf] <- NA
  dAA <- mean(dAA, na.rm = TRUE)
  dBB <- sapply(genesB, function(x) igraph::distances(gr,v=x,to=genesB))
  dBB[dBB == Inf] <- NA
  dBB <- mean(dBB, na.rm = TRUE)
  dAB <- sapply(genesA, function(x) igraph::distances(gr,v=x,to=genesB))
  dAB[dAB == Inf] <- NA
  dAB <- mean(dAB, na.rm = TRUE)
  sAB <- dAB - (dAA + dBB)/2

  return(sAB)
}
