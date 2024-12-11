############################################
## Load Required Libraries
############################################
library(DESeq2)

############################################
## Load Gene Sets
############################################
# Load custom NPC meta-programs
meta_programs <- readxl::read_xlsx("/import/home/lyuah/R/database/meta_programs.xlsx", sheet = 1)
meta_program_list <- as.list(meta_programs)
names(meta_program_list) <- colnames(meta_programs)

# Load EMT gene sets
emt_genes <- readxl::read_xlsx("/import/home/lyuah/R/database/PEMT.xlsx", sheet = 1)
emt_gene_list <- as.list(emt_genes)
names(emt_gene_list) <- colnames(emt_genes)
emt_gene_list <- lapply(emt_gene_list, function(x) x[!is.na(x)])

# Load molecular program gene sets
mp_genes <- read.csv("/import/home/lyuah/R/database/MP_list.csv")
mp_gene_list <- as.list(mp_genes[,-1])

# Load MSigDB gene sets
msigdb_sets <- rbind(
  msigdbr::msigdbr(species="Homo sapiens", category="C2", subcategory = "CP:KEGG"),
  msigdbr::msigdbr(species="Homo sapiens", category="C2", subcategory = "CP:REACTOME"),
  msigdbr::msigdbr(species="Homo sapiens", category="C5", subcategory = "GO:BP")
)
# Process and combine all gene sets
gene_sets <- split(as.character(msigdb_sets$human_gene_symbol), msigdb_sets$gs_name)
gene_sets <- sapply(gene_sets, unique)
gene_sets <- gene_sets[(sapply(gene_sets, length) <= 500) & (sapply(gene_sets, length) > 15)]
gene_sets <- do.call(c, list(mp_gene_list, emt_gene_list, meta_program_list, gene_sets))
names(gene_sets)
# Replace underscores with hyphens in gene set names
names(gene_sets) <- gsub("_", "-", names(gene_sets))

############################################
## Process MCR2017 Dataset
############################################
# Load raw data
mcr2017_counts <- read.delim("/import/home/lyuah/R/NTRK2_PublicData/bulkdata/GSE102349_raw_counts_GRCh38.p13_NCBI.tsv")
gene_annotations <- read.delim("/import/home/lyuah/R/NTRK2_PublicData/bulkdata/Human.GRCh38.p13.annot.tsv")

# Process count matrix
mcr2017_matrix <- as.matrix(mcr2017_counts[,-1])
rownames(mcr2017_matrix) <- gene_annotations$Symbol[match(mcr2017_counts$GeneID, gene_annotations$GeneID)]
mcr2017_filtered <- mcr2017_matrix[rowSums(mcr2017_matrix) >= 50,]

# Normalize with DESeq2
sample_metadata <- data.frame(row.names = colnames(mcr2017_filtered), condition = rep("sample", ncol(mcr2017_filtered)))
mcr2017_dds <- DESeqDataSetFromMatrix(countData = round(mcr2017_filtered), colData = sample_metadata, design = ~ 1)
mcr2017_normalized <- vst(mcr2017_dds)
mcr2017_expression <- assay(mcr2017_normalized)

# Run GSVA analysis
if (!require("GSVA")) { install.packages("GSVA"); library(GSVA) }
mcr2017_gsva <- gsva(mcr2017_expression, gene_sets, method = "gsva", kcdf = "Gaussian", min.sz = 15, max.sz = 500, mx.diff = TRUE, verbose = TRUE)

############################################
## MCR2017 Correlation Analysis
############################################
# Get NTRK2 expression
mcr2017_ntrk2 <- mcr2017_expression["NTRK2",]

# Calculate correlations for each gene
mcr2017_gene_cor <- apply(mcr2017_expression, 1, function(x) {
  cor_test <- cor.test(mcr2017_ntrk2, x, method="pearson")
  return(c(correlation=cor_test$estimate, pvalue=cor_test$p.value))
})

# Create results dataframe for gene-to-gene
mcr2017_gene_cor_results <- data.frame(
  gene = rownames(mcr2017_expression),
  correlation = unlist(mcr2017_gene_cor[1,]),
  pvalue = unlist(mcr2017_gene_cor[2,])
)
mcr2017_gene_cor_results$FDR <- p.adjust(mcr2017_gene_cor_results$pvalue, method="BH")
mcr2017_gene_cor_results <- mcr2017_gene_cor_results[order(abs(mcr2017_gene_cor_results$correlation), decreasing=TRUE),]
head(mcr2017_gene_cor_results)

# Calculate correlations for each pathway
mcr2017_gsva_cor <- apply(mcr2017_gsva, 1, function(x) {
  cor_test <- cor.test(mcr2017_ntrk2, x, method="pearson")
  return(c(correlation=cor_test$estimate, pvalue=cor_test$p.value))
})

# Create results dataframe for pathways
mcr2017_gsva_cor_results <- data.frame(
  geneset = rownames(mcr2017_gsva),
  correlation = unlist(mcr2017_gsva_cor[1,]),
  pvalue = unlist(mcr2017_gsva_cor[2,])
)
mcr2017_gsva_cor_results$FDR <- p.adjust(mcr2017_gsva_cor_results$pvalue, method="BH")
mcr2017_gsva_cor_results <- mcr2017_gsva_cor_results[order(abs(mcr2017_gsva_cor_results$correlation), decreasing=TRUE),]

# Filter out pathways containing NTRK2
mcr2017_gsva_cor_results_filtered <- mcr2017_gsva_cor_results[
  !sapply(gene_sets[mcr2017_gsva_cor_results$geneset], function(x) "NTRK2" %in% x),
]
head(mcr2017_gsva_cor_results_filtered)

############################################
## Process SA2022 Dataset 
############################################
# Load expression data
sa2022_counts <- read.delim("/home/lyuah/R/NTRK2_PublicData/bulkdata/sciadv.abh2445_data_s1.txt")
sa2022_matrix <- as.matrix(sa2022_counts[,-1])
colnames(sa2022_matrix) <- sub("^X", "", colnames(sa2022_matrix))
rownames(sa2022_matrix) <- sa2022_counts[,1]

# Convert gene IDs
if (!require("ensembldb")) BiocManager::install("ensembldb")
if (!require("EnsDb.Hsapiens.v86")) BiocManager::install("EnsDb.Hsapiens.v86")
library(ensembldb)
library(EnsDb.Hsapiens.v86)

ensembl_ids <- sub("\\.[0-9]+$", "", rownames(sa2022_matrix))
edb <- EnsDb.Hsapiens.v86
gene_map <- select(edb, 
  keys = ensembl_ids,
  keytype = "GENEID",
  columns = c("GENEID", "SYMBOL")
)

# Update row names with gene symbols
matrix_rownames <- rownames(sa2022_matrix)
mapped_indices <- match(ensembl_ids, gene_map$GENEID)
valid_symbols <- !is.na(mapped_indices) & gene_map$SYMBOL[mapped_indices] != ""
matrix_rownames[valid_symbols] <- gene_map$SYMBOL[mapped_indices[valid_symbols]]
rownames(sa2022_matrix) <- matrix_rownames

# Filter low-expression genes
sa2022_filtered <- sa2022_matrix[rowSums(sa2022_matrix) >= 50,]

# Normalize with DESeq2
sa2022_metadata <- data.frame(
  row.names = colnames(sa2022_filtered),
  condition = rep("sample", ncol(sa2022_filtered))
)
sa2022_dds <- DESeqDataSetFromMatrix(
  countData = round(sa2022_filtered),
  colData = sa2022_metadata,
  design = ~ 1
)
sa2022_normalized <- vst(sa2022_dds)
sa2022_expression <- assay(sa2022_normalized)

# Run GSVA
sa2022_gsva <- gsva(
  sa2022_expression,
  gene_sets,
  method = "gsva",
  kcdf = "Gaussian",
  min.sz = 15,
  max.sz = 500,
  mx.diff = TRUE,
  verbose = TRUE
)

############################################
## Process Sample Annotations
############################################
sample_info <- read.delim("/home/lyuah/R/NTRK2_PublicData/bulkdata/sciadv.abh2445_data_s2.txt")
cell_info <- read.delim("/home/lyuah/R/NTRK2_PublicData/bulkdata/sciadv.abh2445_data_s3.txt")

library(dplyr)
metadata <- left_join(cell_info, sample_info, by = "patient_index") %>%
  dplyr::select(
    sample_ID, patient_index, dissection_index, cell_type, cell_type_epi,
    recurrence_case, parrafin_type, Gender, Race, EBER, Histo_Pattern,
    Stage, Overall_Stage, Chemo, Chemo_Concurrent, Treatment_Additional, block_age
  )
table(metadata$cell_type)

############################################
## SA2022 Correlation Analysis
############################################
# Filter for tumor samples
tumor_samples <- metadata %>%
  dplyr::filter(cell_type %in% c("TME", "TUM", "NAT", "DYS")) %>%
  pull(sample_ID)
sa2022_expression <- sa2022_expression[, tumor_samples]

# Calculate gene-level correlations
sa2022_gene_cor <- apply(sa2022_expression, 1, function(x) {
  cor_test <- cor.test(sa2022_expression["NTRK2",], x, method="pearson")
  return(c(correlation=cor_test$estimate, pvalue=cor_test$p.value))
})

# Create results dataframe for gene-to-gene
sa2022_gene_cor_results <- data.frame(
  gene = rownames(sa2022_expression),
  correlation = unlist(sa2022_gene_cor[1,]),
  pvalue = unlist(sa2022_gene_cor[2,])
)
sa2022_gene_cor_results$FDR <- p.adjust(sa2022_gene_cor_results$pvalue, method="BH")
sa2022_gene_cor_results <- sa2022_gene_cor_results[order(abs(sa2022_gene_cor_results$correlation), decreasing=TRUE),]
head(sa2022_gene_cor_results, 50)

# Calculate pathway-level correlations
sa2022_gsva <- sa2022_gsva[, tumor_samples]
sa2022_gsva_cor <- apply(sa2022_gsva, 1, function(x) {
  cor_test <- cor.test(sa2022_expression["NTRK2",], x, method="pearson")
  return(c(correlation=cor_test$estimate, pvalue=cor_test$p.value))
})

# Create results dataframe for pathways
sa2022_gsva_cor_results <- data.frame(
  geneset = rownames(sa2022_gsva),
  correlation = unlist(sa2022_gsva_cor[1,]),
  pvalue = unlist(sa2022_gsva_cor[2,])
)
sa2022_gsva_cor_results$FDR <- p.adjust(sa2022_gsva_cor_results$pvalue, method="BH")
sa2022_gsva_cor_results <- sa2022_gsva_cor_results[order(abs(sa2022_gsva_cor_results$correlation), decreasing=TRUE),]

# Filter out pathways containing NTRK2
sa2022_gsva_cor_results_filtered <- sa2022_gsva_cor_results[
  !sapply(gene_sets[sa2022_gsva_cor_results$geneset], function(x) "NTRK2" %in% x),
]
head(sa2022_gsva_cor_results_filtered)

############################################
## Visualization
############################################
library(ggplot2)
library(ggpubr)

# Create correlation scatter plot
mp5_gsva_scores <- sa2022_gsva["MP_5", tumor_samples]
correlation_plot <- ggplot(
  data.frame(NTRK2=sa2022_expression["NTRK2",], MP5=mp5_gsva_scores),
  aes(x=NTRK2, y=MP5)
) +
  geom_point(alpha=0.6) +
  geom_smooth(method="lm", color="red", se=TRUE) +
  labs(
    x="NTRK2 Expression",
    y="MP_5 GSVA Score",
    title="SA Dataset: NTRK2 vs MP_5 Correlation"
  ) +
  theme_bw() +
  stat_cor(method="pearson")
print(correlation_plot)

# Create expression boxplot by cell type
expression_plot <- ggplot(
  data.frame(
    NTRK2 = as.numeric(sa2022_expression["NTRK2",]),
    cell_type = metadata$cell_type[match(colnames(sa2022_expression), metadata$sample_ID)]
  ),
  aes(x = cell_type, y = NTRK2)
) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_bw() +
  labs(
    x = "Cell Type",
    y = "NTRK2 Expression",
    title = "NTRK2 Expression Across Cell Types"
  ) +
  stat_compare_means(method = "anova") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
## Load required libraries and set options
options(future.globals.maxSize = 128000 * 1024^2)
options(Seurat.object.assay.version = "v5")
library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(scCustomize)

## Load reference data
# Load gene annotation data
havana <- read.table("/import/home/lyuah/R/database/HGNC_Mart.txt", header = TRUE, sep = "\t", quote = "", fill = TRUE)

# Load gene length data
gene_info <- read.table("/import/home/lyuah/R/database/gene_pos.txt", header = F, sep = "\t", quote = "", fill = TRUE)
colnames(gene_info) <- c("gene", "chr", "start", "end") 
gene_info$length <- gene_info$end - gene_info$start

## Helper function
counts_to_tpm <- function(counts, gene_lengths) {
  common_genes <- intersect(rownames(counts), gene_lengths$gene)
  counts <- counts[common_genes,]
  matched_lengths <- gene_lengths$length[match(rownames(counts), gene_lengths$gene)]
  rpk <- counts / (matched_lengths/1000)
  scaling_factor <- colSums(rpk)/1e6
  tpm <- sweep(rpk, 2, scaling_factor, "/")
  return(tpm)
}

## Load and process scRNA-seq datasets
# Load NPC tumor cells
npc_tumor <- readRDS("/home/lyuah/R/NTRK2_PublicData/all_npc_tumor_cells.rds")
npc_tumor <- UpdateSeuratObject(npc_tumor)
DefaultAssay(npc_tumor) <- "RNA"
npc_tumor[["SCT"]] <- NULL

# Get count matrix and filter genes
npc_counts <- as.matrix(LayerData(npc_tumor, assay = "RNA", layer = "counts"))
npc_counts <- npc_counts[rowSums(npc_counts > 0) >= 30,]

# Load NPC43 data
npc43 <- readRDS("/home/lyuah/R/NTRK2_PublicData/filted_seu.rds")
npc43 <- UpdateSeuratObject(npc43)
npc43$patient <- "npc43"
npc43$intergration <- npc43$Batch
DefaultAssay(npc43) <- "RNA"
npc43[["SCT"]] <- NULL

# Get count matrix and filter genes
npc43_counts <- as.matrix(LayerData(npc43, assay = "RNA", layer = "counts"))
npc43_counts <- npc43_counts[rowSums(npc43_counts > 0) >= 30,]

# Load C6661 data  
c6661 <- readRDS("/home/lyuah/R/c6661_nomito/Data/Final_C6661_integrated.rds")
c6661 <- UpdateSeuratObject(c6661)
c6661$patient <- "c6661"
c6661$intergration <- c6661$Batch
DefaultAssay(c6661) <- "RNA"
c6661[["SCT"]] <- NULL
c6661[["UCell"]] <- NULL

# Get count matrix and filter genes
c6661_counts <- as.matrix(LayerData(c6661, assay = "RNA", layer = "counts"))
c6661_counts <- c6661_counts[rowSums(c6661_counts > 0) >= 30,]

## Find common genes and filter matrices
common_genes <- Reduce(intersect, list(rownames(npc43_counts), rownames(c6661_counts), rownames(npc_counts)) )

# Filter to common genes
npc43_filt <- npc43_counts[common_genes,]
c6661_filt <- c6661_counts[common_genes,]
npc_filt <- npc_counts[common_genes,]

## Process NPC tumor Seurat object
# Load UCell scores
ucell_scores <- readRDS("/home/lyuah/R/NTRK2_PublicData/Data/UCell_scores.rds")
ucell_scores <- t(ucell_scores)

# Process Seurat object
DefaultAssay(npc_tumor) <- "RNA"
npc_tumor[["RNA"]] <- split(npc_tumor[["RNA"]], f = npc_tumor$intergration)
npc_tumor <- SCTransform(npc_tumor)
npc_tumor <- RunPCA(npc_tumor, npcs = 50, verbose = F)

# Integrate and run UMAP
npc_tumor <- IntegrateLayers(object = npc_tumor, method = CCAIntegration, normalization.method = "SCT", new.reduction = "harmony", verbose = T)
npc_tumor <- RunUMAP(npc_tumor, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony", method = "umap-learn", metric = "correlation")

## Calculate correlations
# Convert counts to TPM
npc43_tpm <- log2(counts_to_tpm(npc43_filt, gene_info) + 1)
c6661_tpm <- log2(counts_to_tpm(c6661_filt, gene_info) + 1)
npc_tpm <- log2(counts_to_tpm(npc_filt, gene_info) + 1)
# Combine TPM matrices
combined_tpm <- cbind(c6661_tpm, npc43_tpm, npc_tpm)

combined_metadata <- bind_rows(npc43@meta.data %>% mutate(dataset = "NPC43") %>% tibble::rownames_to_column("cell_id"), c6661@meta.data %>% mutate(dataset = "C6661") %>% tibble::rownames_to_column("cell_id"), npc_tumor@meta.data %>% mutate(dataset = "Patient") %>% tibble::rownames_to_column("cell_id")) %>% select(cell_id, dataset, intergration, patient)
head(combined_metadata)
table(combined_metadata$dataset)
combined_metadata$intergration[combined_metadata$dataset == "NPC43"] <- "NPC43"
combined_metadata$intergration[combined_metadata$dataset == "C6661"] <- "C6661"
combined_metadata$intergration <- factor(combined_metadata$intergration, levels=names(sort(table(combined_metadata$intergration), decreasing=TRUE)))

# Calculate correlations for each integration
integrations <- unique(combined_metadata$intergration)

# Initialize lists to store correlations
geneset_correlations <- list()
gene_correlations <- list()

for(integ in integrations) {
  cells <- colnames(combined_tpm)[combined_metadata$intergration == integ]
  integ_scores <- ucell_scores[, cells]
  ntrk2_expr <- combined_tpm["NTRK2", cells]
  
  # Calculate gene set correlations
  geneset_corr <- apply(integ_scores, 1, function(x) {
    cor_test <- cor.test(ntrk2_expr, x, method="pearson")
    return(c(correlation=cor_test$estimate, pvalue=cor_test$p.value))
  })
  
  # Calculate gene-to-gene correlations
  gene_corr <- apply(combined_tpm[,cells], 1, function(x) {
    cor_test <- cor.test(ntrk2_expr, x, method="pearson") 
    return(c(correlation=cor_test$estimate, pvalue=cor_test$p.value))
  })
  
  # Create results dataframe for gene sets
  geneset_df <- data.frame(
    geneset = rownames(integ_scores),
    correlation = unlist(geneset_corr[1,]),
    pvalue = unlist(geneset_corr[2,])
  )
  geneset_df$FDR <- p.adjust(geneset_df$pvalue, method="BH")
  geneset_df <- geneset_df[order(abs(geneset_df$correlation), decreasing=TRUE),]
  
  # Create results dataframe for gene-to-gene
  gene_df <- data.frame(
    gene = rownames(combined_tpm),
    correlation = unlist(gene_corr[1,]),
    pvalue = unlist(gene_corr[2,])
  )
  gene_df$FDR <- p.adjust(gene_df$pvalue, method="BH")
  gene_df <- gene_df[order(abs(gene_df$correlation), decreasing=TRUE),]
  
  # Store results in respective lists
  geneset_correlations[[integ]] <- geneset_df
  gene_correlations[[integ]] <- gene_df
}

# Filter out genesets containing NTRK2 for each integration
for(integ in integrations) {
  dim(geneset_correlations[[integ]])
  geneset_correlations[[integ]] <- geneset_correlations[[integ]][
    !sapply(gene_sets[geneset_correlations[[integ]]$geneset], function(x) "NTRK2" %in% x),
  ]
  dim(geneset_correlations[[integ]])
}

# Show filtered results
names(geneset_correlations)
head(geneset_correlations[[1]])
############################################### Intergrate all data

############## List all correlation data

###### Geneset
mcr2017_gsva_cor_results_filtered[1:5,] ## bulk data gene set correlation
sa2022_gsva_cor_results_filtered[1:5,] ## bulk data gene set correlation   
geneset_correlations[[1]][1:5,] ## scRNA-seq gene set correlation

###### Gene
mcr2017_gene_cor_results[1:5,] ## bulk data gene correlation
sa2022_gene_cor_results[1:5,] ## bulk data gene correlation
gene_correlations[[1]][1:5,] ## scRNA-seq gene correlation

############## Select correlation genes

## Select postive correlation genes
# intersect bulk data, for gene correlation from scRNA-seq, the gene should be found in as least 50% of the patients
# 1. intersect bulk data
bulk_genes <- intersect(mcr2017_gene_cor_results$gene[mcr2017_gene_cor_results$FDR < 0.01 & mcr2017_gene_cor_results$correlation > 0.5], sa2022_gene_cor_results$gene[sa2022_gene_cor_results$FDR < 0.01 & sa2022_gene_cor_results$correlation > 0.5])
bulk_genes
# 2. intersect scRNA-seq data
# Get genes that are positively correlated in at least 50% of patients
scrna_genes <- lapply(gene_correlations, function(df) {df$gene[df$FDR < 0.0000001 & df$correlation > 0.2]})
gene_counts <- table(unlist(scrna_genes))
scrna_filtered_genes <- names(gene_counts)[gene_counts >= (length(gene_correlations) * 0.3)]
scrna_filtered_genes
# 3. select postive correlation genes   
postive_correlation_genes <- intersect(bulk_genes, scrna_filtered_genes)
postive_correlation_genes

## Select negative correlation genes
# intersect bulk data, for gene correlation from scRNA-seq, the gene should be found in as least 50% of the patients
# 1. intersect bulk data
bulk_genes <- intersect(mcr2017_gene_cor_results$gene[mcr2017_gene_cor_results$FDR < 0.01 & mcr2017_gene_cor_results$correlation < -0.4], sa2022_gene_cor_results$gene[sa2022_gene_cor_results$FDR < 0.01 & sa2022_gene_cor_results$correlation < -0.4])
bulk_genes
# 2. intersect scRNA-seq data
# Get genes that are positively correlated in at least 50% of patients
scrna_genes <- lapply(gene_correlations, function(df) {df$gene[df$FDR < 0.0000001 & df$correlation < -0.1]})
gene_counts <- table(unlist(scrna_genes))
scrna_filtered_genes <- names(gene_counts)[gene_counts >= (length(gene_correlations) * 0.2)]
scrna_filtered_genes
# 3. select negative correlation genes   
negative_correlation_genes <- bulk_genes
negative_correlation_genes

## Selected genes
selected_genes <- c(postive_correlation_genes, negative_correlation_genes)[-1]
selected_genes

############## Select correlation gene set

## Select postive correlation gene set
# intersect bulk data, for gene correlation from scRNA-seq, the gene should be found in as least 50% of the patients
# 1. intersect bulk data
bulk_genesets <- intersect(mcr2017_gsva_cor_results_filtered$geneset[mcr2017_gsva_cor_results_filtered$FDR < 0.01 & mcr2017_gsva_cor_results_filtered$correlation > 0.3], sa2022_gsva_cor_results_filtered$geneset[sa2022_gsva_cor_results_filtered$FDR < 0.01 & sa2022_gsva_cor_results_filtered$correlation > 0.3])
bulk_genesets
# 2. intersect scRNA-seq data
# Get genes that are positively correlated in at least 50% of patients
scrna_genesets <- lapply(geneset_correlations, function(df) {df$geneset[df$FDR < 0.0000001 & df$correlation > 0.3]})
gene_counts <- table(unlist(scrna_genesets))
scrna_filtered_genesets <- names(gene_counts)[gene_counts >= (length(geneset_correlations) * 0.25)]
scrna_filtered_genesets
# 3. select postive correlation genes   
postive_correlation_genesets <- intersect(bulk_genesets, scrna_filtered_genesets)
postive_correlation_genesets <- postive_correlation_genesets[-3]
postive_correlation_genesets
## Select negative correlation gene set
# intersect bulk data, for gene correlation from scRNA-seq, the gene should be found in as least 50% of the patients
# 1. intersect bulk data
bulk_genesets <- intersect(mcr2017_gsva_cor_results_filtered$geneset[mcr2017_gsva_cor_results_filtered$FDR < 0.01 & mcr2017_gsva_cor_results_filtered$correlation < -0.35], sa2022_gsva_cor_results_filtered$geneset[sa2022_gsva_cor_results_filtered$FDR < 0.01 & sa2022_gsva_cor_results_filtered$correlation < -0.35])
bulk_genesets
# 2. intersect scRNA-seq data
# Get genes that are positively correlated in at least 50% of patients
scrna_genesets <- lapply(geneset_correlations, function(df) {df$geneset[df$FDR < 0.0000001 & df$correlation < -0.3]})
gene_counts <- table(unlist(scrna_genesets))
scrna_filtered_genesets <- names(gene_counts)[gene_counts >= (length(geneset_correlations) * 0.1)]
scrna_filtered_genesets
# 3. select postive correlation genes   
negative_correlation_genesets <- intersect(bulk_genesets, scrna_filtered_genesets)
negative_correlation_genesets <- c("GOBP-POSITIVE-REGULATION-OF-GLYCOLYTIC-PROCESS", "GOBP-MACROAUTOPHAGY")

manually_selected_genesets <- c("GOBP-EPIDERMIS-MORPHOGENESIS", "REACTOME-TRANSCRIPTIONAL-REGULATION-OF-PLURIPOTENT-STEM-CELLS","Hypoxia","Stress")
############## Create correlation plot
## create correlation dataframe
cor_df <- data.frame(Feature = c("SOX2",postive_correlation_genes[-1], negative_correlation_genes, postive_correlation_genesets,manually_selected_genesets, negative_correlation_genesets))
cor_df
# Add correlation results from MCR2017 and SA2022 datasets
cor_df$MCR2017_cor <- sapply(cor_df$Feature, function(f) {
  if(f %in% mcr2017_gene_cor_results$gene) {
    mcr2017_gene_cor_results$correlation[mcr2017_gene_cor_results$gene == f]
  } else if(f %in% mcr2017_gsva_cor_results_filtered$geneset) {
    mcr2017_gsva_cor_results_filtered$correlation[mcr2017_gsva_cor_results_filtered$geneset == f]
  } else {
    NA
  }
})

cor_df$SA2022_cor <- sapply(cor_df$Feature, function(f) {
  if(f %in% sa2022_gene_cor_results$gene) {
    sa2022_gene_cor_results$correlation[sa2022_gene_cor_results$gene == f]
  } else if(f %in% sa2022_gsva_cor_results_filtered$geneset) {
    sa2022_gsva_cor_results_filtered$correlation[sa2022_gsva_cor_results_filtered$geneset == f]
  } else {
    NA
  }
})

# Add correlation results from scRNA-seq data
for(i in 1:length(gene_correlations)){
    cor_df[,i+3] <- sapply(cor_df$Feature, function(f) {
        if(f %in% rownames(gene_correlations[[i]])) {
            gene_correlations[[i]]$correlation[rownames(gene_correlations[[i]]) == f]
        } else if(f %in% geneset_correlations[[i]]$geneset) {
            geneset_correlations[[i]]$correlation[geneset_correlations[[i]]$geneset == f]
        } else {
            NA
        }
    })
    colnames(cor_df)[i+3] <- names(gene_correlations)[i]
}

# Convert correlation dataframe to matrix for heatmap
library(ComplexHeatmap)
cor_matrix <- as.matrix(cor_df[,-1]) # Keep all samples, transpose later
rownames(cor_matrix) <- cor_df$Feature

# Define colors for heatmap
col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))

# Create annotation for feature types and correlation direction
feature_type <- ifelse(grepl("^GO|^MP-|^REACTOME-|Stress|Hypoxia", rownames(cor_matrix)), "Pathway", "Gene")
correlation_direction <- ifelse(rownames(cor_matrix) %in% c(postive_correlation_genes, postive_correlation_genesets, manually_selected_genesets), 
                              "Positive", "Negative")
correlation_direction[1] <- "Positive"
correlation_direction[13] <- "Negative"
correlation_direction[14] <- "Negative"

row_anno <- rowAnnotation(
  Feature_Type = feature_type,
  Correlation = correlation_direction,
  col = list(
    Feature_Type = c("Gene" = "#F8766D", "Pathway" = "#00BFC4"),
    Correlation = c("Positive" = "#EE6677", "Negative" = "#4477AA")
  ),
  show_legend = TRUE
)

# Create column annotation for data types
col_types <- c("Patient Bulk", "Patient Bulk", "Line scRNAseq","Line scRNAseq", rep("Patient scRNA-seq", 14)) # Include all samples
column_anno <- HeatmapAnnotation(
  Data_Type = col_types,
  col = list(Data_Type = c(
    "Patient Bulk" = "#66C2A5",
    "Line scRNAseq" = "#FC8D62", 
    "Patient scRNA-seq" = "#8DA0CB"
  )),
  show_legend = TRUE
)

# Create the heatmap
heatmap <- Heatmap(
  cor_matrix,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8, fontfamily = "Helvetica"),
  column_names_gp = gpar(fontsize = 8),
  top_annotation = column_anno,
  left_annotation = row_anno,
  heatmap_legend_param = list(
    title = "Correlation",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1.0", "-0.5", "0", "0.5", "1.0"),
    color_bar = "continuous",
    legend_height = unit(3, "cm"),
    grid_width = unit(0.3, "cm"),
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 8),
    title_position = "topcenter",
    legend_direction = "vertical",
    legend_side = "right"
  ),
  show_heatmap_legend = TRUE
)

# Save heatmap to PNG file
png("correlation_heatmap.png", width = 1600, height = 1200, res = 200)
draw(heatmap, column_title = "Correlation between NTRK2 and other features in different datasets")
dev.off()

pdf("correlation_heatmap.pdf", width = 10, height = 6)
draw(heatmap, column_title = "Correlation between NTRK2 and other features in different datasets")
dev.off()
