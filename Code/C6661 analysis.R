setwd("/import/home/lyuah/R/c6661")
options(future.globals.maxSize = 128000 * 1024^2)
library(Seurat)
library(dplyr)
library(ggplot2)
library(plot1cell)
library(scCustomize)
library(jsonlite)
library(patchwork)
library(ComplexHeatmap)
library(httpgd)
library(UCell)
# Scanpy default colors
scanpy_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5")
n_colors_needed <- 21
distinct_colors <- colorRampPalette(c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"))(n_colors_needed)
## Load data
havana <- read.table("/import/home/lyuah/R/database/HGNC_Mart.txt", 
                     header = TRUE, 
                     sep = "\t", 
                     quote = "", 
                     fill = TRUE)
mt_genes <- havana$Approved.symbol[
    havana$Chromosome %in% c("mitochondria", "MT", "M", "chrM", "chrMT") & 
    havana$Locus.type != "RNA, ribosomal"]
rRNA <- havana$Approved.symbol[havana$Locus.type == "RNA, ribosomal"]
###################### Direct data analysis (v4 flavor)
## load the 10x output
C6661_UT1 <- Read10X(data.dir = "~/NPC/C6661/10xoutput/UT1_result/outs/filtered_feature_bc_matrix")
C6661_UT2 <- Read10X(data.dir = "~/NPC/C6661/10xoutput/UT2_result/outs/filtered_feature_bc_matrix")
C6661_T241 <- Read10X(data.dir = "~/NPC/C6661/10xoutput/T241_result/outs/filtered_feature_bc_matrix")
C6661_T481 <- Read10X(data.dir = "~/NPC/C6661/10xoutput/T481_result/outs/filtered_feature_bc_matrix")
C6661_T482 <- Read10X(data.dir = "~/NPC/C6661/10xoutput/T482_result/outs/filtered_feature_bc_matrix")
colnames(C6661_UT1) <- paste0("UT1_", colnames(C6661_UT1))
colnames(C6661_UT2) <- paste0("UT2_", colnames(C6661_UT2))
colnames(C6661_T241) <- paste0("T241_", colnames(C6661_T241))
colnames(C6661_T481) <- paste0("T481_", colnames(C6661_T481))
colnames(C6661_T482) <- paste0("T482_", colnames(C6661_T482))
# Combine all data into one matrix
combined_data <- cbind(C6661_UT1, C6661_UT2, C6661_T241, C6661_T481, C6661_T482)
combined_data <- as(combined_data, "dgCMatrix")  # Convert to sparse matrix format
# # Write the counts layer to a directory
library(BPCells)
write_matrix_dir(mat = combined_data, dir = './Data/C6661_counts', overwrite = TRUE)
counts.mat <- open_matrix_dir(dir = './Data/C6661_counts')
# Create a dataframe to annotate cells with their batch
batch_annotation <- data.frame(
  Cell = colnames(combined_data),
  Batch = c(rep("UT1", ncol(C6661_UT1)), rep("UT2", ncol(C6661_UT2)), rep("T241", ncol(C6661_T241)), rep("T481", ncol(C6661_T481)), rep("T482", ncol(C6661_T482)))
)
# Create batch annotation by extracting characters before _ in column names
batch_annotation <- data.frame(
  Cell = colnames(counts.mat),
  Batch = gsub("_.*", "", colnames(counts.mat))
)
head(batch_annotation)
rownames(batch_annotation) <- batch_annotation$Cell
# Create a Seurat v5 object
C6661_combined <- CreateSeuratObject(counts = counts.mat, meta.data = batch_annotation)
C6661_combined$log10counts <- log10(C6661_combined$nCount_RNA)
colnames(C6661_combined@meta.data)
C6661_combined$percent_mito <- PercentageFeatureSet(C6661_combined, features = mt_genes[mt_genes %in% rownames(C6661_combined)])
C6661_combined$percent_ribo <- PercentageFeatureSet(C6661_combined, features = rRNA[rRNA %in% rownames(C6661_combined)])
C6661_combined$percent_mito_ribo <- C6661_combined$percent_mito + C6661_combined$percent_ribo
C6661_combined$percent_ebv <- PercentageFeatureSet(object = C6661_combined, features  = rownames(C6661_combined)[58052:58186])
C6661_combined$Batch_Inter <- "UT"
C6661_combined$Batch_Inter[C6661_combined$Batch == c("T241")] <- "T24"
C6661_combined$Batch_Inter[C6661_combined$Batch %in% c("T481","T482")] <- "T48"
C6661_combined$log10GenesPerUMI <- log10(C6661_combined$nFeature_RNA) / log10(C6661_combined$nCount_RNA)
## filter the data
C6661_filtered <- subset(x = C6661_combined, subset= ((log10counts >= 3.5) & (nFeature_RNA >= 1500)) | (percent_ebv > 0.5))
C6661_filtered <- subset(C6661_filtered, features = names(which(Matrix::rowSums(LayerData(C6661_filtered, layer = "counts")) >= 5)))
C6661_filtered <- subset(x = C6661_combined, subset= ((log10counts >= 3.5) & (nFeature_RNA >= 2000) & (log10GenesPerUMI > 0.8) & (percent_mito < 25) & (percent_mito > 2) & (percent_ribo < 5)) | (percent_ebv > 0.5))
C6661_UT <- subset(C6661_filtered, subset = Batch %in% c("UT1","UT2") )
C6661_T24 <- subset(C6661_filtered, subset = Batch %in% c("T241") )
C6661_T48 <- subset(C6661_filtered, subset = Batch %in% c("T481","T482") )
############ UT data analysis ##########################################################################################################
C6661_UT <- SCTransform(C6661_UT, vst.flavor = "v1", variable.features.n = 4000)
C6661_UT <- CellCycleScoring(C6661_UT, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes,layer = "counts")
C6661_UT <- RunPCA(C6661_UT)
ElbowPlot(C6661_UT, ndims = 50)
C6661_UT <- RunUMAP(C6661_UT, reduction = "pca", dims = 1:20, reduction.name = "umap")
C6661_UT <- FindNeighbors(C6661_UT, reduction = "pca", dims = 1:20)
C6661_UT <- FindClusters(C6661_UT, resolution = 0.5, algorithm = 4)
############ T24 data analysis ##########################################################################################################
C6661_T24 <- SCTransform(C6661_T24, vst.flavor = "v1", variable.features.n = 4000)
C6661_T24 <- CellCycleScoring(C6661_T24, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes,layer = "counts")
C6661_T24 <- RunPCA(C6661_T24)
ElbowPlot(C6661_T24, ndims = 50)
C6661_T24 <- RunUMAP(C6661_T24, reduction = "pca", dims = 1:20, reduction.name = "umap")
C6661_T24 <- FindNeighbors(C6661_T24, reduction = "pca", dims = 1:20)
C6661_T24 <- FindClusters(C6661_T24, resolution = 0.7, algorithm = 4)
############ T48 data analysis ##########################################################################################################
C6661_T48 <- SCTransform(C6661_T48, vst.flavor = "v1", variable.features.n = 4000)
C6661_T48 <- CellCycleScoring(C6661_T48, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes,layer = "counts")
C6661_T48 <- RunPCA(C6661_T48)
ElbowPlot(C6661_T48, ndims = 50)
C6661_T48 <- RunUMAP(C6661_T48, reduction = "pca", dims = 1:30, reduction.name = "umap")
C6661_T48 <- FindNeighbors(C6661_T48, reduction = "pca", dims = 1:30)
C6661_T48 <- FindClusters(C6661_T48, resolution = 0.8, algorithm = 4)
################################### Merge Seurat Objects ###################################
## combine sepreated seurat objects
library(harmony)
colnames(C6661_UT@meta.data)
colnames(C6661_T24@meta.data)
colnames(C6661_T48@meta.data)
colnames(C6661_UT@meta.data)[18] <- "SCT_snn_clusters"
colnames(C6661_T24@meta.data)[18] <- "SCT_snn_clusters"
colnames(C6661_T48@meta.data)[18] <- "SCT_snn_clusters"
## Create integrated object
C6661_integrated <- merge(x = C6661_UT,y = list(C6661_T24, C6661_T48))
C6661_integrated$seurat_clusters2 <- paste0(C6661_integrated$Batch_Inter,"-",C6661_integrated$seurat_clusters)
## Extract UMAP coordinates and add to integrated object
C6661_UT_UMAP <- Embeddings(C6661_UT, "umap")
C6661_T24_UMAP <- Embeddings(C6661_T24, "umap")
C6661_T48_UMAP <- Embeddings(C6661_T48, "umap")
C6661_Seperation_UMAP <- do.call(rbind, list(C6661_UT_UMAP, C6661_T24_UMAP, C6661_T48_UMAP))
head(C6661_Seperation_UMAP)
tail(C6661_Seperation_UMAP)
dim(C6661_Seperation_UMAP)
C6661_integrated[['SepUMAP']] <- CreateDimReducObject(embeddings = C6661_Seperation_UMAP, key = "SepUMAP_", global = T, assay = "RNA")
C6661_integrated[['umap']] <- CreateDimReducObject(embeddings = C6661_Seperation_UMAP, key = "UMAP_", global = T, assay = "RNA")
################################### Caculate UCell scores ###################################
pan_nmf <- readxl::read_xlsx("/import/home/lyuah/R/database/meta_programs.xlsx",sheet = 1)
pan_nmf_list <- as.list(pan_nmf)
names(pan_nmf_list) <- colnames(pan_nmf)
emt <- readxl::read_xlsx("/import/home/lyuah/R/database/PEMT.xlsx",sheet = 1)
emt_list <- as.list(emt)
names(emt_list) <- colnames(emt)
emt_list <- lapply(emt_list, function(x) x[!is.na(x)])
npc_list <- read.csv("/import/home/lyuah/R/database/MP_list.csv")
npc_list <- as.list(npc_list[,-1])
###### GO meta program analysis
h.human <- rbind(msigdbr::msigdbr(species="Homo sapiens",category="C2",subcategory = "CP:KEGG"),
                 #msigdbr(species="Homo sapiens",category="C2",subcategory = "CP:PID"),
                 msigdbr::msigdbr(species="Homo sapiens",category="C2",subcategory = "CP:REACTOME"),
                 msigdbr::msigdbr(species="Homo sapiens",category="C5",subcategory = "GO:BP")
                 )
h.sets = split(as.character(h.human$human_gene_symbol), h.human$gs_name)
h.sets = sapply(h.sets, unique)
h.sets = h.sets[(sapply(h.sets, length) <= 500) & (sapply(h.sets, length) >15)]
h.sets = do.call(c, list(npc_list,emt_list,pan_nmf_list, h.sets))
C6661_Ucell_scores <- ScoreSignatures_UCell(as.matrix(LayerData(C6661_preliminary_filtered, assay = "RNA", layer = "counts")), features=h.sets, BPPARAM=BiocParallel::SnowParam(24))
C6661_integrated[["UCell"]] <- CreateAssayObject(data = t(UCell_scores))
DefaultAssay(C6661_integrated) <- "RNA"
C6661_integrated[['SCT']] = NULL
C6661_integrated[['RNA']]$data = NULL
C6661_integrated[['RNA']]$scale.data = NULL
colnames(C6661_integrated@meta.data)
C6661_integrated$Batch_Inter <- factor(C6661_integrated$Batch_Inter, levels = c("UT","T24","T48"))
C6661_integrated <- SCTransform(C6661_integrated, vst.flavor = "v1", variable.features.n = 4000)
C6661_integrated <- RunPCA(C6661_integrated)
## harmony integration
C6661_integrated <- RunHarmony(C6661_integrated, "Batch", nclust = 40, max_iter = 10, early_stop = T)
ElbowPlot(C6661_integrated, ndims = 50,reduction = "harmony")
C6661_integrated <- RunUMAP(C6661_integrated, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", method = "umap-learn", metric = "correlation", min.dist = 0.4, spread = 2)
## annotation
C6661_integrated$States2[C6661_integrated$Batch_Inter == "T48" & C6661_integrated$seurat_clusters == "1"] <- "T48_Cytolysis"
C6661_integrated$States2[C6661_integrated$Batch_Inter == "T48" & C6661_integrated$seurat_clusters == "5"] <- "T48_G2M"
C6661_integrated$States2[C6661_integrated$Batch_Inter == "T48" & C6661_integrated$seurat_clusters == "12"] <- "T48_Lytic"
C6661_integrated$States2[C6661_integrated$Batch_Inter == "T48" & C6661_integrated$seurat_clusters %in% c("9","7","10")] <- "T48_PreLytic"
C6661_integrated$States2[C6661_integrated$Batch_Inter == "T48" & C6661_integrated$seurat_clusters %in% c("2","3","4","6","8","11")] <- "T48_NR_like"

C6661_integrated$States2[C6661_integrated$Batch_Inter == "T24" & C6661_integrated$seurat_clusters %in% c("5","2","3")] <- "T24_Intermediate"
C6661_integrated$States2[C6661_integrated$Batch_Inter == "T24" & C6661_integrated$seurat_clusters %in% c("4")] <- "T24_PreLytic"
C6661_integrated$States2[C6661_integrated$Batch_Inter == "T24" & C6661_integrated$Phase == "G2M"] <- "T24_G2M"
C6661_integrated$States2[C6661_integrated$Batch_Inter == "T24" & C6661_integrated$Phase == "S"] <- "T24_S"

C6661_integrated$States2[C6661_integrated$seurat_clusters %in% c("1","7")] <- "UT_Epi_like"
C6661_integrated$States2[C6661_integrated$seurat_clusters %in% c("5","2","8") & (C6661_integrated$Phase != "G2M")] <- "UT_NR_like"
C6661_integrated$States2[C6661_integrated$Phase == "G2M"] <- "UT_G2M"
C6661_integrated$States2[C6661_integrated$Phase == "S"] <- "UT_S"
## DE; similar for GO and genes
C6661_UT_temp <- subset(C6661_integrated, Batch_Inter == "UT")
C6661_T24_temp <- subset(C6661_integrated, Batch_Inter == "T24")
C6661_T48_temp <- subset(C6661_integrated, Batch_Inter == "T48")
table(Idents(C6661_UT_temp))
table(Idents(C6661_T24_temp))
table(Idents(C6661_T48_temp))
C6661_UT_DE <- FindAllMarkers(C6661_UT_temp, min.pct = 0.1, logfc.threshold = 0.25)
C6661_T24_DE <- FindAllMarkers(C6661_T24_temp, min.pct = 0.1, logfc.threshold = 0.25)
C6661_T48_DE <- FindAllMarkers(C6661_T48_temp, min.pct = 0.1, logfc.threshold = 0.25)
DefaultAssay(C6661_UT_Ucell) <- "UCell"
DefaultAssay(C6661_T24_Ucell) <- "UCell"
DefaultAssay(C6661_T48_Ucell) <- "UCell"
UT_GO_Ucell <- FindAllMarkers(C6661_UT_Ucell)
T24_GO_Ucell <- FindAllMarkers(C6661_T24_Ucell)
T48_GO_Ucell <- FindAllMarkers(C6661_T48_Ucell)
## also caculate the specificity score for GO
UT_GO_RSS <- SCENIC::calcRSS(AUC=LayerData(C6661_UT_Ucell, assay = "UCell"), cellAnnotation=C6661_UT_Ucell$seurat_clusters)
UT_GO_RSS <- as.data.frame(as.table(UT_GO_RSS), responseName = "AUcell") %>% group_by(Var2)
colnames(UT_GO_RSS) <- c("GO","Cluster","RSS")
############# gene analysis
markers <- C6661_UT_DE %>% dplyr::filter(cluster == "NR_like") %>% 
  mutate(min_nonzero_p = min(p_val_adj[p_val_adj > 0], na.rm = TRUE),
         log_p = case_when(p_val_adj == 0 ~ -log10(min_nonzero_p),TRUE ~ -log10(p_val_adj)),
         combined_score = (log_p)^(1/3) * avg_log2FC * abs(pct.1 - pct.2)) %>% arrange(desc(combined_score))
head(markers, n = 30)
FeaturePlot_scCustom(C6661_UT_temp, features = markers$gene[1:6], colors_use = c('#4477AA','lemonchiffon1','indianred1','darkred'), reduction = "SepUMAP",num_columns = 3)
FeaturePlot_scCustom(C6661_T24_temp, features = markers$gene[1:6], colors_use = c('#4477AA','lemonchiffon1','indianred1','darkred'), reduction = "SepUMAP",num_columns = 3)
FeaturePlot_scCustom(C6661_T48_temp, features = markers$gene[1:6], colors_use = c('#4477AA','lemonchiffon1','indianred1','darkred'), reduction = "SepUMAP",num_columns = 3)
############# pathway analysis
####### G2M
UT_G2M_GO <- inner_join(UT_GO_Ucell %>% dplyr::filter(cluster == "G2M"), UT_GO_RSS %>% dplyr::filter(cluster == "G2M"), by = c("gene","cluster")) %>% 
  mutate(min_nonzero_p = min(p_val_adj[p_val_adj > 0], na.rm = TRUE),
         log_p = case_when(p_val_adj == 0 ~ -log10(min_nonzero_p),TRUE ~ -log10(p_val_adj)),
         combined_score = (log_p)^(1/3) * avg_log2FC * RSS * (pct.1 - pct.2)) %>% arrange(desc(combined_score))
head(UT_G2M_GO)
T24_G2M_GO <- inner_join(T24_GO_Ucell %>% dplyr::filter(cluster == "G2M"), T24_GO_RSS %>% dplyr::filter(cluster == "G2M"), by = c("gene","cluster")) %>% 
  mutate(min_nonzero_p = min(p_val_adj[p_val_adj > 0], na.rm = TRUE),
         log_p = case_when(p_val_adj == 0 ~ -log10(min_nonzero_p),TRUE ~ -log10(p_val_adj)),
         combined_score = (log_p)^(1/3) * avg_log2FC * RSS * (pct.1 - pct.2)) %>% arrange(desc(combined_score))
head(T24_G2M_GO)
T48_G2M_GO <- inner_join(T48_GO_Ucell %>% dplyr::filter(cluster == "G2M"), T48_GO_RSS %>% dplyr::filter(cluster == "G2M"), by = c("gene","cluster")) %>% 
  mutate(min_nonzero_p = min(p_val_adj[p_val_adj > 0], na.rm = TRUE),
         log_p = case_when(p_val_adj == 0 ~ -log10(min_nonzero_p),TRUE ~ -log10(p_val_adj)),
         combined_score = (log_p)^(1/3) * avg_log2FC * RSS * (pct.1 - pct.2)) %>% arrange(desc(combined_score))
head(T48_G2M_GO)
G2M_GO <- Reduce(intersect, list(UT_G2M_GO$gene[1:10],T24_G2M_GO$gene[1:10],T48_G2M_GO$gene[1:10]))
G2M_GO
########################################## Plot cluster to cluster Correlation
cor_features = Reduce(intersect, list(C6661_UT@assays$SCT@var.features[1:2650], C6661_T24@assays$SCT@var.features[1:2650], C6661_T48@assays$SCT@var.features[1:2650]))
cor_seurat = merge(x = C6661_UT,y = list(C6661_T24, C6661_T48))
print(length(cor_features))
cor_seurat$States_anno = paste0(C6661_integrated$Batch_Inter,"_",C6661_integrated$cellular_states2)
pseudo_bulk = as.matrix(AggregateExpression(cor_seurat, group.by = "States_anno", assays = "SCT", slot = "scale.data", verbose = T, features = cor_features)$SCT)

pdf(paste0("./Revision_plot/Revision_SI_003_correlation.pdf"), width=6, height=4)
col <- colorRampPalette(c("#4477AA","#77AADD","#FFFFFF","#EE9988", "#BB4444"))
corrplot::corrplot(cor(pseudo_bulk, use = "complete.obs", method = "spearman"), order = 'hclust', hclust.method = "average", method = 'color', col = col(200),  addCoef.col = "black", tl.col = "black", tl.srt = 45, number.cex = 0.48, tl.cex = 0.8, rect.lwd = 0.5, cl.ratio = 0.25, cl.cex = 0.75)
dev.off()          
