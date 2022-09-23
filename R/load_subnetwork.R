#' Loads the subnetwork information from pairs of drugs and diseases
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
#' load_subnetwork('Myeloid Leukemia, Chronic', 'Multiple Sclerosis', 'Imatinib', 'Tandutinib')
#' }

load_subnetwork <- function(DISEASE1, DISEASE2, DRUG1, DRUG2){
  utils::data("disgenet")
  utils::data("interactome_subset")
  # Filter disgenet df by score
  disgenet <- disgenet %>%
    filter(score >= 0.3)

  # Convert interactome to graph
  SUBNETWORK <- igraph::graph_from_data_frame(interactome_subset)
  SUBNETWORK <- igraph::as.undirected(SUBNETWORK) #convert to undirected
  SUBNETWORK <- igraph::simplify(SUBNETWORK, remove.multiple = TRUE, remove.loops = TRUE)

  # Find genes and subset disgenet and interactome
  disease1_genes <- as.character(disgenet[disgenet$diseaseName %in% DISEASE1,]$geneSymbol)
  disease2_genes <- as.character(disgenet[disgenet$diseaseName %in% DISEASE2,]$geneSymbol)
  drug1_genes <- as.character(interactome_subset[interactome_subset$Protein_A == DRUG1,]$Protein_B)
  drug2_genes <- as.character(interactome_subset[interactome_subset$Protein_A == DRUG2,]$Protein_B)
  disease1_genes <- disease1_genes[which(disease1_genes %in% igraph::V(SUBNETWORK)$name)]
  disease2_genes <- disease2_genes[which(disease2_genes %in% igraph::V(SUBNETWORK)$name)]

  SUBNETWORK <- rbind(cbind(interactome_subset[(interactome_subset$Protein_A%in% disease1_genes),], source = DISEASE1),
                  cbind(interactome_subset[(interactome_subset$Protein_B%in% disease1_genes),], source = DISEASE1),
                  cbind(interactome_subset[(interactome_subset$Protein_A%in% disease2_genes),], source = DISEASE2),
                  cbind(interactome_subset[(interactome_subset$Protein_B%in% disease2_genes),], source = DISEASE2),
                  cbind(interactome_subset[interactome_subset$Protein_A %in% DRUG1,], source = DRUG1),
                  cbind(interactome_subset[interactome_subset$Protein_A %in% DRUG2,], source = DRUG2))
  # Remove duplicate entries from subnetwork df
  SUBNETWORK <- SUBNETWORK[!duplicated(SUBNETWORK[,c("Protein_A", "Protein_B")]),]

  SUBNETGRAPH <- igraph::graph_from_data_frame(SUBNETWORK)
  SUBNETGRAPH <- igraph::as.undirected(SUBNETGRAPH) #convert to undirected
  SUBNETGRAPH <- igraph::simplify(SUBNETGRAPH, remove.multiple = TRUE, remove.loops = TRUE)

  # Gene sources
  geneSources <- as.data.frame(rbind(cbind(genes = disease1_genes, source = 1),
                      cbind(genes = disease2_genes, source = 2),
                      cbind(genes = drug1_genes, source = 3),
                      cbind(genes = drug2_genes, source = 4)))
  geneSources.wide <- tidyr::pivot_wider(geneSources, names_from = "source", values_from = "source")
  geneSources.wide$source <- apply(geneSources.wide[ ,c("1","2","3","4")], 1, function(x) paste(x[!is.na(x)], collapse = ","))
  geneSources.wide <- geneSources.wide[,c("genes", "source")]

  # Add source of gene as node attributes
  SUBNETGRAPH <- set_vertex_attr(SUBNETGRAPH, "source", index = geneSources.wide$genes, geneSources.wide$source)
  # SUBNETGRAPH <- delete.vertices(SUBNETGRAPH, V(SUBNETGRAPH)[which(is.na(get.vertex.attribute(SUBNETGRAPH, "source")))])

  return(list(subnetwork=SUBNETGRAPH,
              disease1_genes=disease1_genes, disease2_genes=disease2_genes,
              drug1_genes=drug1_genes, drug2_genes=drug2_genes))
}

# igraph::betweenness(SUBNETGRAPH, v=DRUG1, directed=FALSE)
# igraph::betweenness(SUBNETGRAPH, v=DRUG2, directed=FALSE)
