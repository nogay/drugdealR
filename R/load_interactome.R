#' Loads the interactome information from data
#'
#' @param DISEASE1 name of the first disease. Should be compatible with
#' DisNor (https://disnor.uniroma2.it/) disease names.
#' @param DISEASE2 name of the second disease. Should be compatible with
#' DisNor (https://disnor.uniroma2.it/) disease names.
#' @param DRUG1 name of the first drug. Should be compatible with Drugbank
#' (https://go.drugbank.com/) drug names.
#' @param DRUG2 name of the second drug. Should be compatible with Drugbank
#' (https://go.drugbank.com/) drug names.
#'
#' @return list
#' @export
#'
#' @examples
#' \dontrun{
#' library(drugdealR)
#' load_interactome('Myeloid Leukemia, Chronic', 'Multiple Sclerosis', 'Imatinib', 'Tandutinib')
#' }

load_interactome <- function(DISEASE1, DISEASE2, DRUG1, DRUG2){
  utils::data("disgenet")
  utils::data("interactome_subset")

  SUBNETWORK <- igraph::graph_from_data_frame(interactome_subset)
  SUBNETWORK <- igraph::as.undirected(SUBNETWORK) #convert to undirected
  SUBNETWORK <- igraph::simplify(SUBNETWORK, remove.multiple = TRUE, remove.loops = TRUE)

  disease1_genes <- as.character(disgenet[disgenet$diseaseName %in% DISEASE1,]$geneSymbol)
  disease2_genes <- as.character(disgenet[disgenet$diseaseName %in% DISEASE2,]$geneSymbol)
  drug1_genes <- as.character(interactome_subset[interactome_subset$Protein_A == DRUG1,]$Protein_B)
  drug2_genes <- as.character(interactome_subset[interactome_subset$Protein_A == DRUG2,]$Protein_B)
  disease1_genes <- disease1_genes[which(disease1_genes %in% igraph::V(SUBNETWORK)$name)]
  disease2_genes <- disease2_genes[which(disease2_genes %in% igraph::V(SUBNETWORK)$name)]

  SUBNETWORK <- rbind(interactome_subset[(interactome_subset$Protein_A%in% disease1_genes),],
                  interactome_subset[(interactome_subset$Protein_B%in% disease1_genes),],
                  interactome_subset[(interactome_subset$Protein_A%in% disease2_genes),],
                  interactome_subset[(interactome_subset$Protein_B%in% disease2_genes),],
                  interactome_subset[interactome_subset$Protein_A %in% DRUG1,],
                  interactome_subset[interactome_subset$Protein_A %in% DRUG2,])

  INTERACTOME <- igraph::graph_from_data_frame(SUBNETWORK)
  INTERACTOME <- igraph::as.undirected(INTERACTOME)

  SUBNETGRAPH <- igraph::graph_from_data_frame(SUBNETWORK)
  SUBNETGRAPH <- igraph::as.undirected(SUBNETGRAPH) #convert to undirected
  SUBNETGRAPH <- igraph::simplify(SUBNETGRAPH, remove.multiple = TRUE, remove.loops = TRUE)

  return(list(interactome=INTERACTOME, subnetwork=SUBNETGRAPH,
              disease1_genes=disease1_genes, disease2_genes=disease2_genes,
              drug1_genes=drug1_genes, drug2_genes=drug2_genes))
}

#measure_separation(disease1_genes, disease2_genes, INTERACTOME)
#measure_separation(drug1_genes, drug2_genes, INTERACTOME)

# igraph::betweenness(SUBNETGRAPH, v=DRUG1, directed=FALSE)
# igraph::betweenness(SUBNETGRAPH, v=DRUG2, directed=FALSE)
