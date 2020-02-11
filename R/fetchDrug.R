#' Load Drugbank data (downloaded January 2020)
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
  DRUG_DF$genes <- gsub("\\'", "", DRUG_DF$genes)

  DRUG_DF <- separate(DRUG_DF, 'genes', paste0('genes', c(1:(max(str_count(DRUG_DF$genes, ","))+1))), sep = ',', remove = T)

  DRUG_DF <- pivot_longer(DRUG_DF,
    cols = starts_with("genes"),
    names_to = "genes",
    values_drop_na = TRUE
  )
  DRUG_DF <- DRUG_DF[DRUG_DF$value != "",]
  DRUG_DF <- DRUG_DF[,c("name","approval","value")]

  DRUG_DF <- separate(DRUG_DF, "value", c("gene","action"), sep=":")

  return(DRUG_DF)
}
