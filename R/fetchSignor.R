#' Load SIGNOR 2.0 data (Downloaded February 2020)
#'
#' @return dataframe
#' @export
#'
#' @examples
#' fetchSignor()
fetchSignor <- function(){
  SIGNOR_DF <- read.csv("./data/signor_11022020.tsv", sep="\t")
  SIGNOR_DF <- SIGNOR_DF[,c("ENTITYA","TYPEA","IDA","ENTITYB","TYPEB","IDB","EFFECT")]

  return(SIGNOR_DF)
}
