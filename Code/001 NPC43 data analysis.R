## Load library
library(Seurat)
library(scCustomize)
library(ggsci)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(singleseqgset)
library(UCell)
##### Filter QC
npc43_seu <- CreateSeuratObject(npc43_mat, meta.data = npc43_anno)
npc43_seu$log10GenesPerUMI <- log10(npc43_seu$nFeature_RNA) / log10(npc43_seu$nCount_RNA)
npc43_seu$percent.mt <- PercentageFeatureSet(object = npc43_seu, pattern = "^MT-")
VlnPlot(npc43_seu, features = c("nCount_RNA","nFeature_RNA","percent.mt", "log10GenesPerUMI"), group.by = "Batch", pt.size = 0.1, ncol = 4)
filtered_npc <- subset(x = npc43_seu, subset= (nCount_RNA >= 1000) & (nFeature_RNA >= 1000) & (log10GenesPerUMI > 0.75) & (log10GenesPerUMI < 0.87) & (percent.mt < 25) & (percent.mt > 4) )
filtered_npc <- subset(filtered_npc, features = rownames(filtered_npc)[Matrix::rowSums(filtered_npc) > 5] )
dim(filtered_npc);dim(npc43_seu)
VlnPlot(filtered_npc, features = c("nCount_RNA","nFeature_RNA","percent.mt", "log10GenesPerUMI"), group.by = "Batch", pt.size = 0.1, ncol = 4)
filtered_npc <- NormalizeData(filtered_npc)
filtered_npc <- ScaleData(filtered_npc)
filtered_npc <- CellCycleScoring(filtered_npc,s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
##### UMAP
filtered_npc <- SCTransform(filtered_npc, vst.flavor = "v1", variable.features.n = 4000)
filtered_npc <- RunPCA(filtered_npc, npcs = 50, verbose = FALSE)
ElbowPlot(filtered_npc,ndims = 50)
filtered_npc <- RunUMAP(filtered_npc,dims = c(1:13), seed.use = 888)

DimPlot(filtered_npc, group.by = "Batch", pt.size = 1) + scale_color_locuszoom()
DimPlot(filtered_npc, group.by = "Phase", pt.size = 1) + scale_color_locuszoom()
FeaturePlot(filtered_npc, "NTRK2", pt.size = 1) + scale_colour_gradientn(colours =c("#023e8a", "#e9ecef", "#e63946"))
FeaturePlot(filtered_npc, "log10GenesPerUMI", pt.size = 1) + scale_colour_gradientn(colours =c("#023e8a", "#e9ecef", "#e63946"))
FeaturePlot(filtered_npc, c("nCount_RNA","nFeature_RNA","percent.mt", "log10GenesPerUMI"), pt.size = 1)
##### Non-responsive cells
filtered_npc = FindNeighbors(filtered_npc, reduction = "umap", dims = 1:2, verbose = FALSE)
filtered_npc = FindClusters(filtered_npc, resolution = 0.01, verbose = FALSE, algorithm =1)
DimPlot(filtered_npc)

filtered_npc$Batch2 <- "npc43"
filtered_npc$Batch2[filtered_npc$Batch == "npct0"] <- "Untreated NPC43"
filtered_npc$Batch2[filtered_npc$Batch == "npct0-2"] <- "Untreated NPC43"
filtered_npc$Batch2[filtered_npc$Batch == "npct24"] <- "Treated NPC43 (24hr)"
filtered_npc$Batch2[filtered_npc$Batch == "npct48-1"] <- "Treated NPC43 (48hr)"
filtered_npc$Batch2[filtered_npc$Batch == "npct48-2"] <- "Treated NPC43 (48hr)"
DimPlot(filtered_npc,group.by = "Batch2", pt.size = 0.1) + mytheme2 + ggtitle(label = NULL) + scale_color_manual(values = c("#0072B2","#3B4992FF", "#EEA236FF" ))

filtered_npc$nonresponsive <- "Induced state"
nonresponsive_cells <- colnames(filtered_npc)[filtered_npc$seurat_clusters == 1][colnames(filtered_npc)[filtered_npc$seurat_clusters == 1] %in% colnames(filtered_npc)[filtered_npc$Batch2 %in% c("Treated NPC43 (24hr)","Treated NPC43 (48hr)")]]
native_cells = colnames(filtered_npc)[filtered_npc$seurat_clusters == 1][colnames(filtered_npc)[filtered_npc$seurat_clusters == 1] %in% colnames(filtered_npc)[filtered_npc$Batch2 %in% c("Untreated NPC43")]]
filtered_npc$nonresponsive[colnames(filtered_npc) %in% nonresponsive_cells] <- "Nonresponsive state"
filtered_npc$nonresponsive[colnames(filtered_npc) %in% native_cells] <- "Native state"
DimPlot(filtered_npc,group.by = "nonresponsive", pt.size = 0.1) + scale_color_manual(values = c( "#3B4992FF","#EEA236FF", "#BB0021FF")) + mytheme2 + ggtitle(label = NULL)

library(ComplexHeatmap)
pdf(paste0("S-001-nonresponsive_heatmap_",1000, "_wardD_euclidean.pdf"), width = 10, height = 6)
plot_hvg_heatmap(filtered_npc, feature_num = 1000, clustering_method = "ward.D", clustering_dist = "euclidean")
dev.off()
################################### Cellular states identification and annotation
## UT cells clustering
# 1. dimension reduction
UT_npc <- subset(filtered_npc, subset = (Batch %in% c("npct0-2", "npct0") ) )
DefaultAssay(UT_npc) <- "RNA"
UT_npc[['SCT']] = NULL
UT_npc[['RNA']]$data = NULL
UT_npc[['RNA']]$scale.data = NULL
UT_npc[["RNA"]] <- split(UT_npc[["RNA"]], f = UT_npc$Batch)
UT_npc <- NormalizeData(UT_npc)
UT_npc <- ScaleData(UT_npc)
UT_npc <- SCTransform(UT_npc, vst.flavor = "v1", variable.features.n = 4000, vars.to.regress = c("log10GenesPerUMI"))
UT_npc <- RunPCA(UT_npc, npcs = 100, verbose = FALSE)
UT_npc <- RunUMAP(UT_npc, dims = c(1:30), seed.use = 66, spread = 2)
# 2. clustering
UT_npc = FindNeighbors(UT_npc, reduction = "pca", dims = 1:30, verbose = TRUE)
UT_npc = FindClusters(UT_npc, resolution = 0.5, verbose = TRUE, algorithm =4)
UT_npc$cellular_states <- "Others"
UT_npc$cellular_states[UT_npc$Phase == "S"] <- "Cycling-S"
UT_npc$cellular_states[UT_npc$Phase == "G2M"] <- "Cycling-G2M"
UT_npc$cellular_states[UT_npc$Phase == "G1" & UT_npc$SCT_snn_res.0.3 == "2"] <- "UT-Cluster-1"
UT_npc$cellular_states[UT_npc$Phase == "G1" & UT_npc$SCT_snn_res.0.3 == "1"] <- "UT-Cluster-2"
UT_npc$cellular_states[UT_npc$Phase == "G1" & UT_npc$SCT_snn_res.0.3 %in% c("3","4")] <- "UT-Cluster-3"
DimPlot(UT_npc,group.by = "cellular_states", pt.size = 0.1) + scale_color_lancet()+ mytheme2 + ggtitle(label = NULL)
# 3. DEGs
UT_DEs <- FindAllMarkers(UT_npc_nocc, only.pos = TRUE, min.pct = 0.1)

## T24 cells clustering
# 1. dimension reduction
T24_npc <- subset(filtered_npc, subset = (Batch %in% c("npct24") ) )
DefaultAssay(T24_npc) <- "RNA"
T24_npc[['SCT']] = NULL
T24_npc[['RNA']]$data = NULL
T24_npc[['RNA']]$scale.data = NULL
T24_npc <- NormalizeData(T24_npc)
T24_npc <- ScaleData(T24_npc)
T24_npc <- SCTransform(T24_npc, vst.flavor = "v1", variable.features.n = 4000)
T24_npc <- RunPCA(T24_npc, npcs = 50, verbose = FALSE)
T24_npc <- RunUMAP(T24_npc,dims = c(1:20), seed.use = 66, min.dist = 0.3)
# 2. clustering
T24_npc = FindNeighbors(T24_npc, reduction = "pca", dims = 1:20, verbose = FALSE)
T24_npc = FindClusters(T24_npc, resolution = 0.5, verbose = FALSE, algorithm =4)
T24_npc$cellular_states <- "Others"
T24_npc$cellular_states[T24_npc$Phase == "S"] <- "Cycling-S"
T24_npc$cellular_states[T24_npc$Phase == "G2M"] <- "Cycling-G2M"
T24_npc$cellular_states[T24_npc$Phase == "G1" & T24_npc$SCT_snn_res.0.5 == "8"] <- "T24-Cluster-1"
T24_npc$cellular_states[T24_npc$Phase == "G1" & T24_npc$SCT_snn_res.0.5 == "1"] <- "T24-Cluster-3"
T24_npc$cellular_states[T24_npc$Phase == "G1" & T24_npc$SCT_snn_res.0.5 %in% c("3","5")] <- "T24-Cluster-2"
T24_npc$cellular_states <- factor(T24_npc$cellular_states, levels = c("T24-Cluster-1","T24-Cluster-2","T24-Cluster-3","Cycling-S","Cycling-G2M","Others"))
DimPlot(T24_npc,group.by = "cellular_states", pt.size = 0.1) + scale_color_lancet()+ mytheme2 + ggtitle(label = NULL)
# 3. DEGs
T24_DEs <- FindAllMarkers(UT_npc_nocc, only.pos = TRUE, min.pct = 0.1)

## T48 cells clustering
# 1. dimension reduction
T48_npc <- subset(npc43_seu_v5, subset = (Batch %in% c("npct48-1", "npct48-2") ) )
DefaultAssay(T48_npc) <- "RNA"
T48_npc[['SCT']] = NULL
T48_npc[['RNA']]$data = NULL
T48_npc[['RNA']]$scale.data = NULL
T48_npc <- NormalizeData(T48_npc)
T48_npc <- ScaleData(T48_npc)
T48_npc <- SCTransform(T48_npc, vst.flavor = "v1", variable.features.n = 4000, vars.to.regress = c("percent.mt", "log10GenesPerUMI"))
T48_npc <- RunPCA(T48_npc, npcs = 100, verbose = FALSE)
T48_npc <- RunUMAP(T48_npc,dims = c(1:20), seed.use = 8)
# 2. clustering
T48_npc = FindNeighbors(T48_npc, reduction = "pca", dims = 1:16, verbose = FALSE)
T48_npc = FindClusters(T48_npc, resolution = 0.4, verbose = FALSE, algorithm =4)
T48_npc$cellular_states <- "Others"
T48_npc$cellular_states[T48_npc$Phase == "S"] <- "Cycling-S"
T48_npc$cellular_states[T48_npc$Phase == "G2M"] <- "Cycling-G2M"
T48_npc$cellular_states[T48_npc$Phase == "G1" & T48_npc$SCT_snn_res.0.35 %in% c("4", "1")] <- "T48-Cluster-3"
T48_npc$cellular_states[T48_npc$Phase == "G1" & T48_npc$SCT_snn_res.0.35 %in% c("5", "3")] <- "T48-Cluster-2"
T48_npc$cellular_states[T48_npc$SCT_snn_res.0.35 == "7"] <- "T48-Cluster-1"
T48_npc$cellular_states[T48_npc$SCT_snn_res.0.35 == "8"] <- "T48-Cluster-4"
T48_npc$cellular_states <- factor(T48_npc$cellular_states, levels = c("T48-Cluster-1","T48-Cluster-2","T48-Cluster-3","Cycling-S","Cycling-G2M","Others","T48-Cluster-4"))
DimPlot(T48_npc,group.by = "cellular_states", pt.size = 0.1) + scale_color_lancet()+ mytheme2 + ggtitle(label = NULL)
# 3. DEGs
T48_DEs <- FindAllMarkers(UT_npc_nocc, only.pos = TRUE, min.pct = 0.1)

########### intergration with seurat and harmony
npc43_seu_v5 = filtered_npc
npc43_seu_v5[['SCT']] = NULL
npc43_seu_v5[['RNA']]$data = NULL
npc43_seu_v5[['RNA']]$scale.data = NULL
npc43_seu_v5[["RNA"]] <- split(npc43_seu_v5[["RNA"]], f = npc43_seu_v5$Batch)
npc43_seu_v5 <- NormalizeData(npc43_seu_v5)
npc43_seu_v5 <- ScaleData(npc43_seu_v5)
npc43_seu_v5 <- SCTransform(npc43_seu_v5, vst.flavor = "v1", variable.features.n = 4000, vars.to.regress = c("log10GenesPerUMI"))
npc43_seu_v5 <- RunPCA(npc43_seu_v5, npcs = 100, verbose = FALSE)
npc43_seu_v5 <- CellCycleScoring(npc43_seu_v5,s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
npc43_seu_v5 <- RunUMAP(npc43_seu_v5, dims = c(1:20), seed.use = 666, min.dist = 0.4, spread = 2)

npc43_seu_v5 <- IntegrateLayers(npc43_seu_v5,  method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", assay = "SCT", dims = 1:20)
npc43_seu_v5 <- RunUMAP(npc43_seu_v5, reduction = "harmony", dims = 1:20, return.model = T, verbose = T, min.dist = 0.2)
########### Score cells with different genesets
h.human <- rbind(#msigdbr::msigdbr(species="Homo sapiens",category="C2",subcategory = "CP:KEGG"),
                 #msigdbr(species="Homo sapiens",category="C2",subcategory = "CP:PID"),
                 #msigdbr::msigdbr(species="Homo sapiens",category="C2",subcategory = "CP:REACTOME"),
                 msigdbr::msigdbr(species="Homo sapiens",category="C5",subcategory = "GO:BP")
                 )
h.sets = split(as.character(h.human$human_gene_symbol), h.human$gs_name)
h.sets = sapply(h.sets, unique)
h.sets = h.sets[(sapply(h.sets, length) <= 500) & (sapply(h.sets, length) >15)]
## can download cancer meta program from: "https://www.dropbox.com/scl/fi/22xtcdh0z7bnn5g5ugz33/meta_programs_2023-07-13.xlsx?rlkey=2e7d718s46zybiyvjuptm67n4&dl=1"
npc_ucell_scores <- ScoreSignatures_UCell(GetAssayData(npc43_seu_v5, assay = "SCT", slot = "data"), features=h.sets, name = "")
npc43_seu_v5[["Ucell"]] = CreateAssayObject(t(npc_ucell_scores))
DefaultAssay(npc43_seu_v5) = "Ucell"
npc43_seu_v5_Ucell = FindAllMarkers(npc43_seu_v5, assay = "Ucell", logfc.threshold = 0, test.use = "wilcox", only.pos = TRUE)
########### Score cells with PROGENy
library(progeny)
CellsClusters <- data.frame(Cell = names(Idents(npc43_seu_v5)),
                            CellType = paste0(npc43_seu_v5$Batch3),
                            stringsAsFactors = FALSE)

pre_progeny <- progeny(npc43_seu_v5, scale=FALSE, organism="Human", verbose = T,assay_name = "RNA",
                       top=100, perm=1,return_assay = TRUE)

pre_progeny <- Seurat::ScaleData(pre_progeny, assay = "progeny")
progeny_scores_df <-
  as.data.frame(t(GetAssayData(pre_progeny, slot = "scale.data",
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

progeny_scores_df$CellType = factor(progeny_scores_df$CellType, levels= c("UT","T24","T48"))

summarized_progeny_scores <- progeny_scores_df %>%
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
paletteLength <- 100
myColor <- colorRampPalette(c("Darkblue", "white","red"))(paletteLength)
progenyBreaks <- c(seq(min(summarized_progeny_scores_df), 0,
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(summarized_progeny_scores_df)/paletteLength,
                       max(summarized_progeny_scores_df),
                       length.out=floor(paletteLength/2)))

pheatmap(t(summarized_progeny_scores_df),fontsize=10,
         fontsize_row = 10,cluster_cols = F, cluster_rows = F,
         color=myColor, breaks = progenyBreaks,
         main = "10X PROGENy (500)", angle_col = 90,
         treeheight_col = 0,  border_color = NA)

stat.test <- progeny_scores_df %>% filter(Pathway %in% c("EGFR", "MAPK","PI3K","TGFb")) %>%
  group_by(Pathway) %>%
  rstatix::t_test(Activity ~ CellType, ref.group = "UT") %>%
  rstatix::adjust_pvalue(method = "bonferroni") %>%
  rstatix::add_significance()
stat.test <- stat.test %>% rstatix::add_xy_position(x = "CellType")

ggplot(aes(x= CellType, y = Activity, color = Pathway, fill = Pathway), data = progeny_scores_df %>% filter(Pathway %in% c("EGFR", "MAPK","PI3K","TGFb")) ) +
  geom_violin(position=position_dodge(),draw_quantiles=c(0.5)) +
  geom_boxplot(width=0.2, lwd = 0.1, color="black",position = position_dodge(width =0.9), outlier.size = 0.001) +
  scale_fill_brewer(palette="RdBu") +   scale_color_brewer(palette="RdBu") +
  facet_wrap(vars(`Pathway`),nrow = 1) +
  theme_bw() + theme(axis.text = element_text(colour = "black", family = "sans",size = 12)) +
  stat_pvalue_manual(stat.test)
ggsave(filename = "PROGENy selected pathway scores across Treatment.svg",device = svg, units = "cm",width = 16, height = 6)
########### Some example code to produce figure in paper
## Plot individual UMAP
plot_1 = DimPlot(UT_npc,group.by = "native_cellular_states", pt.size = 0.1) +
  scale_color_manual(values = color_pal[c(1,2,3,6,7,8)] ) +
  mytheme_single_umap + ggtitle(label = "Untreated NPC43")
plot_2 = DimPlot(T24_npc,group.by = "native_cellular_states", pt.size = 0.1) +
  scale_color_manual(values = color_pal[c(1,2,4,5,6,8)] ) +
  mytheme_single_umap + ggtitle(label = "Treated NPC43 (24hr)")
plot_3 = DimPlot(T48_npc,group.by = "native_cellular_states", pt.size = 0.1) +
  scale_color_manual(values = color_pal[c(1,2,4,5,6,7,8,9)] ) +
  mytheme_single_umap + ggtitle(label = "Treated NPC43 (48hr)")
ggpubr::ggarrange(plot_1,plot_2,plot_3, ncol = 3,legend = "right",widths = c(1,1.1,1.1))
## Plot correlation heatmap
cor_features = Reduce(intersect, list(UT_npc@assays$SCT@var.features[1:320],
                                      T24_npc@assays$SCT@var.features[1:320],
                                      T48_npc@assays$SCT@var.features[1:320]
                                      )
                     )
pseudo_bulk = as.matrix(AggregateExpression(npc43_seu_v5,group.by = "test",assays = "SCT",slot = "scale.data",verbose = T, features = cor_features)$SCT)
corrplot::corrplot(cor(pseudo_bulk,use = "complete.obs",method = "spearman"),
         col=col(200), order = 'hclust', hclust.method = "average", method = 'color',addrect = 3,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, number.cex = 0.7)
## plot selected markers with dot plot
complex_dotplot_multiple(npc43_seu_v5, features = c("SOX2","NTRK2","BCL2","DST","CAV1","MMP10","HSPA5","S100A9"),
                         group = "Batch3",
                         celltypes =  c("State-01","State-1.5","State-02","State-03","State-04"))
## volcano plot and surface markers
volcano_Data = T48_DEs %>% dplyr::filter(cluster == "State-01(NTRK2-high)")
volcano_Data$pct.dif = volcano_Data$pct.1 - volcano_Data$pct.2
volcano_Data$log10_pval = -log10(volcano_Data$p_val_adj)
volcano_Data$Fold_changes = abs(volcano_Data$avg_log2FC)
ggplot(volcano_Data) + geom_point(aes(x = pct.1, y = pct.dif, color = log10_pval,size = Fold_changes)) +
  geom_text_repel(data = volcano_Data %>%
                    dplyr::filter(gene %in% c("NTRK2", "SOX2", "BCL2", "ADAM23", "FOXE1","GAST")),
                  aes(label = gene, x = pct.1, y = pct.dif), box.padding = unit(0.3, "lines"),hjust= 0.30, min.segment.length	= 0,size = 3.5) &
  scale_colour_gradientn(colours =c("#023e8a", "lemonchiffon1", "darkred")) & theme_cowplot(font_family = "Arial")
