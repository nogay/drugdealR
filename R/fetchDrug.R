#' Load Drugbank data (January 2020)
#'
#' @return dataframe
#' @export
#'
#' @examples
#' fetchDrug()

fetchDrugbank <- function(){
  DRUG_DF <- read.csv(paste0("./data/drugbank.tsv"), sep = "\t", stringsAsFactors = FALSE)
  colnames(DRUG_DF) <- c("name", "approval", "genes")
  DRUG_DF$genes <- gsub("\\{", "", DRUG_DF$genes)
  DRUG_DF$genes <- gsub("\\}", "", DRUG_DF$genes)

  DRUG_DF <- DRUG_DF %>% separate('genes', paste0('genes', c(1:(max(str_count(df$genes, ","))+1))), sep = ',', remove = T)

  return(DRUG_DF)
}
