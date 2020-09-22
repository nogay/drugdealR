#' Loads the interactome information from data
#'
#' @return dataframe
#' @export
#'
#' @examples
#' int_df <- load_interactome()
#'

load_interactome <- function(){
  DISEASE1 <- "Myeloid Leukemia, Chronic"
  DISEASE2 <- "Multiple Sclerosis"
  DRUG1 <- "Imatinib"
  DRUG2 <- "Tandutinib"

  # Load interactome and make graph ###################
  interactome_subset <- read.csv("./data/megazord.csv")
  interactome_subset$X <- NULL
  disgenet <- read.csv("./data/disgenet.tsv", sep="\t")
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
}



#measure_separation(disease1_genes, disease2_genes, INTERACTOME)
#measure_separation(drug1_genes, drug2_genes, INTERACTOME)

# igraph::betweenness(SUBNETGRAPH, v=DRUG1, directed=FALSE)
# igraph::betweenness(SUBNETGRAPH, v=DRUG2, directed=FALSE)
