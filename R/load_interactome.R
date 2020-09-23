#' Loads the interactome information from data
#'
#' @param DISEASE1
#' @param DISEASE2
#' @param DRUG1
#' @param DRUG2
#'
#' @return list
#' @export
#'
#' @examples
#' library(drugdealR)
#' load_interactome()

load_interactome <- function(DISEASE1, DISEASE2, DRUG1, DRUG2){
  interactome_subset <- utils::read.csv("./data/megazord.csv")
  interactome_subset$X <- NULL
  disgenet <- utils::read.csv("./data/disgenet.tsv", sep="\t")
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

  return(list(interactome=INTERACTOME, subnetwork=SUBNETGRAPH))
}

#measure_separation(disease1_genes, disease2_genes, INTERACTOME)
#measure_separation(drug1_genes, drug2_genes, INTERACTOME)

# igraph::betweenness(SUBNETGRAPH, v=DRUG1, directed=FALSE)
# igraph::betweenness(SUBNETGRAPH, v=DRUG2, directed=FALSE)
