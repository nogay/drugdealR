#' Load Drugbank data
#'
#' @param x character
#'
#' @return A dataframe
#' @export
#'
#' @examples
#' fetchDrug("./data/drugbank.xml")
fetchDrugbank <- function(x){
  xmldoc <- xmlParse(x)
  rootNode <- xmlRoot(xmldoc)
  data <- xmlSApply(rootNode,function(x) xmlSApply(x, xmlValue))
  cd.catalog <- data.frame(t(data),row.names=NULL)
}
