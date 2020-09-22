measure_dists <- function(disease, drug, gr, subnet){
  gr_noDrug <- subnet[-which(subnet$Protein_A == drug),]
  gr_noDrug <- igraph::graph_from_data_frame(gr_noDrug)
  gr_noDrug <- igraph::as.undirected(gr_noDrug) #convert to undirected
  gr_noDrug <- igraph::simplify(gr_noDrug, remove.multiple = TRUE, remove.loops = TRUE)
  disNetDists <- sapply(disease, function(x) igraph::distances(gr,v=x,to=disease))
  disNetDists <- as.numeric(disNetDists)
  disNetDists_noDrug <- sapply(disease, function(x) igraph::distances(gr_noDrug,v=x,to=disease))
  rownames(disNetDists_noDrug) <- disease
  disNetDists_noDrug <- as.numeric(disNetDists_noDrug)
  return(all.equal(disNetDists,disNetDists_noDrug))
}

measure_separation <- function(diseaseA, diseaseB, gr){
  dAA <- sapply(diseaseA, function(x) igraph::distances(gr,v=x,to=diseaseA))
  dAA[dAA == Inf] <- NA
  dAA <- mean(dAA, na.rm = TRUE)
  dBB <- sapply(diseaseB, function(x) igraph::distances(gr,v=x,to=diseaseB))
  dBB[dBB == Inf] <- NA
  dBB <- mean(dBB, na.rm = TRUE)
  dAB <- sapply(diseaseA, function(x) igraph::distances(gr,v=x,to=diseaseB))
  dAB[dAB == Inf] <- NA
  dAB <- mean(dAB, na.rm = TRUE)
  sAB <- dAB - (dAA + dBB)/2

  return(sAB)
}
