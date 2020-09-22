library(igraph)

# Load interactome and make graph ###################
interactome_subset <- read.csv("./data/megazord.csv")
interactome_subset$X <- NULL
gr <- graph_from_data_frame(interactome_subset)
gr <- as.undirected(gr) #convert to undirected
gr <- simplify(gr, remove.multiple = TRUE, remove.loops = TRUE)

# Get shortest path distances single disease ##############
imatinib <- interactome_subset[interactome_subset$Protein_A == "Imatinib",]
natalizumab <- interactome_subset[interactome_subset$Protein_A == "Natalizumab",]
tandutinib <- interactome_subset[interactome_subset$Protein_A == "Tandutinib",]
MYELOIDLEUKEMIA <- c("HSPA5","GSTP1","BIRC5","ALOX5","ABL1","PRAME","SETBP1","CA9","PER1","BCR","PER3","PER2","IFNA2","ARNTL","RUNX1","HMMR","MAPK14","PRTN3","CRY1","CRY2","WT1")
MYELOIDLEUKEMIA <- MYELOIDLEUKEMIA[which(MYELOIDLEUKEMIA %in% V(gr)$name)]
MULTIPLESCLEROSIS <- c("AHI1","SPEF2", "KIF1B", "CLDN11", "SLC30A7","HLA-DRA", "IL7", "CBLB", "P2RX7", "TYK2", "TAGAP","VAV2", "CD24", "HLA-DQB1", "PVR", "MERTK", "ZFP36L1","CD86", "EPS15L1", "SP140", "CNR1", "VCAM1", "TNFSF14","INAVA", "MAPK1", "ZMIZ1", "HACE1", "IFNB1", "POMC","RPS6KB1", "CD6", "MPHOSPH9", "MCAM", "CIITA", "ARHGEF10","ZNF433", "HLA-F", "STAT3", "MLANA", "BCHE", "IL7R", "MGAT5","SELE", "CLSTN2", "TNFRSF1A", "CLEC16A", "VSIG8", "PDCD1", "BTNL2", "APOE", "SLC15A2", "EVI5", "ODF3B", "TIMMDC1","MALT1", "IRF8", "CYP24A1", "IL2RA", "MPV17L2", "PRKRA","MANBA", "CXCR5", "BATF", "CYP27B1", "IL1B", "METTL1", "BACH2","RAD21L1", "ZBTB46", "SAE1", "GPC5", "CD58", "TSBP1", "AGAP2","SUMF1", "DLEU1", "SLC11A1", "NCOA5", "HLA-DPB1", "DKKL1","KCNJ10", "IL12A", "CLECL1", "CD40", "HLA-DRB1", "RNASEL","ICAM1", "SYK", "STAT4", "MYNN", "PTPRC", "ZNF767P", "RGS14", "TNFAIP3")
MULTIPLESCLEROSIS <- MULTIPLESCLEROSIS[which(MULTIPLESCLEROSIS %in% V(gr)$name)]

SUBNET <- rbind(interactome_subset[(interactome_subset$Protein_A%in% MULTIPLESCLEROSIS),],
                interactome_subset[(interactome_subset$Protein_B%in% MULTIPLESCLEROSIS),],
                interactome_subset[(interactome_subset$Protein_A%in% MYELOIDLEUKEMIA),],
                interactome_subset[(interactome_subset$Protein_B%in% MYELOIDLEUKEMIA),],
                interactome_subset[interactome_subset$Protein_A == "Imatinib",],
                interactome_subset[interactome_subset$Protein_A == "Natalizumab",],
                interactome_subset[interactome_subset$Protein_A == "Tandutinib",])
gr <- graph_from_data_frame(SUBNET)
gr <- as.undirected(gr) #convert to undirected
gr <- simplify(gr, remove.multiple = TRUE, remove.loops = TRUE)
plot(gr, #edge.color=ecol,
     edge.width=2, #edge.arrow.mode=0,
     vertex.size = 5, edge.arrow.size=.2,
     vertex.label = NA)

# no Imatinib Sclerosis
gr_noImatinib <- SUBNET[-which(SUBNET$Protein_A == "Imatinib"),]
gr_noImatinib <- graph_from_data_frame(gr_noImatinib)
gr_noImatinib <- as.undirected(gr_noImatinib) #convert to undirected
gr_noImatinib <- simplify(gr_noImatinib, remove.multiple = TRUE, remove.loops = TRUE)
SCLEROSISDISTS <- sapply(MULTIPLESCLEROSIS, function(x) distances(gr,v=x,to=MULTIPLESCLEROSIS))
SCLEROSISDISTS <- as.numeric(SCLEROSISDISTS)
SCLEROSISDISTS_noImatinib <- sapply(MULTIPLESCLEROSIS, function(x) distances(gr_noImatinib,v=x,to=MULTIPLESCLEROSIS))
rownames(SCLEROSISDISTS_noImatinib) <- MULTIPLESCLEROSIS
SCLEROSISDISTS_noImatinib <- as.numeric(SCLEROSISDISTS_noImatinib)
all.equal(SCLEROSISDISTS,SCLEROSISDISTS_noImatinib)

# no Imatinib Myeloidleukemia
gr_noImatinib <- SUBNET[-which(SUBNET$Protein_A == "Imatinib"),]
gr_noImatinib <- graph_from_data_frame(gr_noImatinib)
gr_noImatinib <- as.undirected(gr_noImatinib) #convert to undirected
gr_noImatinib <- simplify(gr_noImatinib, remove.multiple = TRUE, remove.loops = TRUE)
MYELOIDLEUKEMIADISTS <- sapply(MYELOIDLEUKEMIA, function(x) distances(gr,v=x,to=MYELOIDLEUKEMIA))
MYELOIDLEUKEMIADISTS <- as.numeric(MYELOIDLEUKEMIADISTS)
MYELOIDLEUKEMIADISTS_noImatinib <- sapply(MYELOIDLEUKEMIA, function(x) distances(gr_noImatinib,v=x,to=MYELOIDLEUKEMIA))
rownames(MYELOIDLEUKEMIADISTS_noImatinib) <- MYELOIDLEUKEMIA
MYELOIDLEUKEMIADISTS_noImatinib <- as.numeric(MYELOIDLEUKEMIADISTS_noImatinib)
all.equal(MYELOIDLEUKEMIADISTS,MYELOIDLEUKEMIADISTS_noImatinib)

# no Natalizumab Sclerosis
gr_noNatalizumab <- SUBNET[-which(SUBNET$Protein_A == "Natalizumab"),]
gr_noNatalizumab <- graph_from_data_frame(gr_noNatalizumab)
gr_noNatalizumab <- as.undirected(gr_noNatalizumab) #convert to undirected
gr_noNatalizumab <- simplify(gr_noNatalizumab, remove.multiple = TRUE, remove.loops = TRUE)
SCLEROSISDISTS <- sapply(MULTIPLESCLEROSIS, function(x) distances(gr,v=x,to=MULTIPLESCLEROSIS))
SCLEROSISDISTS <- as.numeric(SCLEROSISDISTS)
SCLEROSISDISTS_noNatalizumab <- sapply(MULTIPLESCLEROSIS, function(x) distances(gr_noNatalizumab,v=x,to=MULTIPLESCLEROSIS))
rownames(SCLEROSISDISTS_noNatalizumab) <- MULTIPLESCLEROSIS
SCLEROSISDISTS_noNatalizumab <- as.numeric(SCLEROSISDISTS_noNatalizumab)
all.equal(SCLEROSISDISTS,SCLEROSISDISTS_noNatalizumab)

# no Natalizumab Myeloidleukemia
gr_noNatalizumab <- SUBNET[-which(SUBNET$Protein_A == "Natalizumab"),]
gr_noNatalizumab <- graph_from_data_frame(gr_noNatalizumab)
gr_noNatalizumab <- as.undirected(gr_noNatalizumab) #convert to undirected
gr_noNatalizumab <- simplify(gr_noNatalizumab, remove.multiple = TRUE, remove.loops = TRUE)
MYELOIDLEUKEMIADISTS <- sapply(MYELOIDLEUKEMIA, function(x) distances(gr,v=x,to=MYELOIDLEUKEMIA))
MYELOIDLEUKEMIADISTS <- as.numeric(MYELOIDLEUKEMIADISTS)
MYELOIDLEUKEMIADISTS_noNatalizumab <- sapply(MYELOIDLEUKEMIA, function(x) distances(gr_noNatalizumab,v=x,to=MYELOIDLEUKEMIA))
rownames(MYELOIDLEUKEMIADISTS_noNatalizumab) <- MYELOIDLEUKEMIA
MYELOIDLEUKEMIADISTS_noNatalizumab <- as.numeric(MYELOIDLEUKEMIADISTS_noNatalizumab)
all.equal(MYELOIDLEUKEMIADISTS,MYELOIDLEUKEMIADISTS_noNatalizumab)

# no Tandutinib Sclerosis
gr_noTandutinib <- SUBNET[-which(SUBNET$Protein_A == "Tandutinib"),]
gr_noTandutinib <- graph_from_data_frame(gr_noTandutinib)
gr_noTandutinib <- as.undirected(gr_noTandutinib) #convert to undirected
gr_noTandutinib <- simplify(gr_noTandutinib, remove.multiple = TRUE, remove.loops = TRUE)
SCLEROSISDISTS <- sapply(MULTIPLESCLEROSIS, function(x) distances(gr,v=x,to=MULTIPLESCLEROSIS))
SCLEROSISDISTS <- as.numeric(SCLEROSISDISTS)
SCLEROSISDISTS_noTandutinib <- sapply(MULTIPLESCLEROSIS, function(x) distances(gr_noTandutinib,v=x,to=MULTIPLESCLEROSIS))
rownames(SCLEROSISDISTS_noTandutinib) <- MULTIPLESCLEROSIS
SCLEROSISDISTS_noTandutinib <- as.numeric(SCLEROSISDISTS_noTandutinib)
all.equal(SCLEROSISDISTS,SCLEROSISDISTS_noTandutinib)

# no Tandutinib Myeloidleukemia
gr_noTandutinib <- SUBNET[-which(SUBNET$Protein_A == "Tandutinib"),]
gr_noTandutinib <- graph_from_data_frame(gr_noTandutinib)
gr_noTandutinib <- as.undirected(gr_noTandutinib) #convert to undirected
gr_noTandutinib <- simplify(gr_noTandutinib, remove.multiple = TRUE, remove.loops = TRUE)
MYELOIDLEUKEMIADISTS <- sapply(MYELOIDLEUKEMIA, function(x) distances(gr,v=x,to=MYELOIDLEUKEMIA))
MYELOIDLEUKEMIADISTS <- as.numeric(MYELOIDLEUKEMIADISTS)
MYELOIDLEUKEMIADISTS_noTandutinib <- sapply(MYELOIDLEUKEMIA, function(x) distances(gr_noTandutinib,v=x,to=MYELOIDLEUKEMIA))
rownames(MYELOIDLEUKEMIADISTS_noTandutinib) <- MYELOIDLEUKEMIA
MYELOIDLEUKEMIADISTS_noTandutinib <- as.numeric(MYELOIDLEUKEMIADISTS_noTandutinib)
all.equal(MYELOIDLEUKEMIADISTS,MYELOIDLEUKEMIADISTS_noTandutinib)
