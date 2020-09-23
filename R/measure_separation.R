#' Measure separation of two entities (disease or drugs)
#'
#' @param entityA name of the first item to measure separation. Should be compatible
#' with DisNor (https://disnor.uniroma2.it/) disease names. If is a drug,
#' should be compatible with Drugbank names.
#' @param entityB name of the second item to measure separation. Should be compatible
#' with DisNor (https://disnor.uniroma2.it/) disease names if is a disease. If is a drug,
#' should be compatible with Drugbank names.
#' @param gr the subnetwork igraph object that contains all genes involved in
#' diseases (according to DisNor) and drug targets (according to Drugbank) of
#' interest.
#'
#' @return data.frame
#' @export
#'
#' @examples
#' \dontrun{
#' library(drugdealR)
#' load_interactome()
#'
#' To calculate separation between two drugs, insert drug names
#' measure_separation('Imatinib', 'Tandutinib', gr)
#'
#' To calculate separation between two diseases, insert disease names
#' measure_separation('Myeloid Leukemia, Chronic', 'Multiple Sclerosis')
#' }

measure_separation <- function(entityA, entityB, gr){
  dAA <- sapply(entityA, function(x) igraph::distances(gr,v=x,to=entityA))
  dAA[dAA == Inf] <- NA
  dAA <- mean(dAA, na.rm = TRUE)
  dBB <- sapply(entityB, function(x) igraph::distances(gr,v=x,to=entityB))
  dBB[dBB == Inf] <- NA
  dBB <- mean(dBB, na.rm = TRUE)
  dAB <- sapply(entityA, function(x) igraph::distances(gr,v=x,to=entityB))
  dAB[dAB == Inf] <- NA
  dAB <- mean(dAB, na.rm = TRUE)
  sAB <- dAB - (dAA + dBB)/2

  return(sAB)
}
