## After finish the SENIC in python, we can import to R to finish the downstream analysis
# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
# For some of the plots:
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
# load loom file in R
scenicLoomPath1 <- file.path("/import/home/lyuah/Py_content/npc_patient_senic/npc43_scenic_out.loom")
loom1 <- open_loom(scenicLoomPath1)
    # Read information from loom file:
    regulons_incidMat1 <- get_regulons(loom1, column.attr.name="Regulons")
    regulons1 <- regulonsToGeneLists(regulons_incidMat1)
    regulonAUC1 <- get_regulons_AUC(loom1, column.attr.name="RegulonsAUC")
    regulonAucThresholds1 <- get_regulon_thresholds(loom1)
close_loom(loom1)
write.csv(getAUC(regulonAUC2),"npc_scenic_genie3_auc.csv")
reg_data1 = read.table("/import/home/lyuah/Py_content/npc_patient_senic/grnboost2_adj.tsv")
############ SOX2 TARGET NTRK2 analysis
ntrk2_tfs_genie3 = names(regulons2)[which(sapply(regulons2, function(x) "NTRK2" %in% x))]
#################### use network matrix to plot network
## 1. TF selection; select TF specificly activated in state-01 cells
## use genie3 result
cellClusters = data.frame(row.names = colnames(npc43_seu_v5), type = npc43_seu_v5$test)
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,"type"])
regulonAUC <- regulonAUC2[onlyNonDuplicatedExtended(rownames(regulonAUC2)),]
regulonActivity_byCellType <- sapply(cellsPerCluster, function(cells) rowMeans(getAUC(regulonAUC2)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled = regulonActivity_byCellType_Scaled[,]
hm = draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity", row_names_gp=grid::gpar(fontsize=6))); hm # row font size
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
rss <- calcRSS(AUC=getAUC(regulonAUC2), cellAnnotation=cellClusters[colnames(regulonAUC2), selectedResolution])
plotRSS_oneSet2(rss, setName = "UT State-01", n = 10)
########## extract the top TFs for NR cells; only check TFs activated in NR cells
selected_tf = unique(c(getTOP_tfs(rss, "UT State-01", n = 25), getTOP_tfs(rss, "T24 State-01", n = 20), getTOP_tfs(rss, "T48 State-01", n = 20)))
selected_tf = substr(selected_tf, 1, nchar(selected_tf) - 3)
regulons_matrix  = as.matrix(regulons_incidMat2)
rownames(regulons_matrix) = substr(rownames(regulons_incidMat2), 1, nchar(rownames(regulons_incidMat2)) - 3)
regulons_matrix = regulons_matrix[,c(rownames(regulons_matrix), "NTRK2")]
regulons_matrix = rbind(regulons_matrix, rep(0,dim(regulons_matrix)[2]))
rownames(regulons_matrix)[dim(regulons_matrix)[2]] = "NTRK2"
g = regulons_matrix[c("NTRK2",selected_tf),c("NTRK2",selected_tf)]
g  <- graph_from_adjacency_matrix(g, mode = "directed")
paths <- all_simple_paths(g, "SOX2", "NTRK2") # This step may take long time if the adjacency_matrix is large
nodes <- unique(unlist(paths))
sub_net <- subgraph(g, nodes)
sub_net <- simplify(sub_net, remove.multiple = FALSE, remove.loops = TRUE)
write.csv(as.matrix(sub_net),"net_matrix.csv", quote=FALSE) ## This martrix can be import to Cytoscape for plotting and analysis
