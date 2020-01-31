#' Load Drugbank data (January 2020)
#'
#' @return dataframe
#' @export
#'
#' @examples
#' fetchDrug()

fetchDrugbank <- function(){
  DRUG_DF <- read.csv(paste0("./data/drugbank.tsv"), sep = "\t")
}
