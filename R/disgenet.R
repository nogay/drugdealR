#' DisGeNet data (https://www.disgenet.org/). Contains disease related genes.
#'
#' Includes various properties of evidence for a given gene-disease pair.
#'
#' @format A data frame with 84038 rows and 16 variables:
#' \describe{
#'  \item{geneId}{ID of gene}
#'  \item{geneSymbol}{Scientific symbol of gene}
#'  \item{DSI}{Disease Specificity Index for the gene}
#'  \item{DPI}{Disease Pleitropy Index for the gene}
#'  \item{diseaseId}{ID of disease}
#'  \item{diseaseName}{Name of disease}
#'  \item{diseaseType}{Disease type. One of: phenotype, group, disease}
#'  \item{diseaseClass}{Class of disease}
#'  \item{diseaseSemanticType}{Detailed type}
#'  \item{score}{Gene-Disease Association score}
#'  \item{EI}{Evidence Index}
#'  \item{YearInitial}{Year of first evidence}
#'  \item{YearFinal}{Year of last evidence}
#'  \item{NofPmids}{Number of PMIDs supporting association}
#'  \item{NofSnps}{Number of SNPs for this Gene-Disease association}
#'  \item{source}{Source of evidence}
#' }
"disgenet"
