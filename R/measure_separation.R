#' Measure separation of two entities (disease or drugs)
#'
#' @param diseaseA
#' @param diseaseB
#' @param gr
#'
#' @return data.frame
#' @export
#'
#' @examples
#' library(drugdealR)
#' load_interactome()
#'
#' To calculate separation between two drugs, insert drug names
#' measure_separation('Imatinib', 'Tandutinib', gr)
#'
#' To calculate separation between two diseases, insert disease names
#' measure_separation('Myeloid Leukemia, Chronic', 'Multiple Sclerosis')

measure_separation <- function(diseaseA, diseaseB, gr){
  dAA <- sapply(diseaseA, function(x) igraph::distances(gr,v=x,to=diseaseA))
  dAA[dAA == Inf] <- NA
  dAA <- mean(dAA, na.rm = TRUE)
  dBB <- sapply(diseaseB, function(x) igraph::distances(gr,v=x,to=diseaseB))
  dBB[dBB == Inf] <- NA
  dBB <- mean(dBB, na.rm = TRUE)
  dAB <- sapply(diseaseA, function(x) igraph::distances(gr,v=x,to=diseaseB))
  dAB[dAB == Inf] <- NA
  dAB <- mean(dAB, na.rm = TRUE)
  sAB <- dAB - (dAA + dBB)/2

  return(sAB)
}
