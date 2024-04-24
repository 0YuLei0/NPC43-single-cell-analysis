## color_pal
color_pal = c( "#E7AB9A","#EEA236FF","#587850", "#D80032", "#FF6969","#9632B8FF", "#6499E9", "#12486B" ,"#BB2525")
## plot themes
mytheme2 <- theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
                  panel.background = element_rect(fill = "white", colour = "white"),
                  panel.grid = element_blank(),
                  axis.title = element_text(size=8, family = "Arial", colour = "black", face = "bold"),
                  axis.text.x = element_text(size=8, family = "Arial", colour = "black", face = "plain"),
                  axis.text.y = element_text(size=8, family = "Arial", colour = "black", face = "plain"),
                  legend.title = element_blank(),
                  legend.text = element_text(size=8, family = "Arial", face = "plain", colour = "black"),
                  legend.key = element_rect(fill = NA)
)
mytheme <- theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
                 panel.background = element_rect(fill = "white", colour = "white"),
                 panel.grid = element_blank(),
                 axis.title = element_text(size=12, family = "Arial", colour = "black", face = "bold"),
                 axis.text.x = element_text(size=8, family = "Arial", colour = "black", face = "plain", angle = 30,hjust = 1),
                 axis.text.y = element_text(size=10, family = "Arial", colour = "black", face = "plain"),
                 legend.title = element_blank(),
                 legend.text = element_text(size=8, family = "Arial", face = "plain", colour = "black"),
                 legend.key = element_rect(fill = NA)
)
## Plot Heatmap with hihgly vriable genes
plot_hvg_heatmap <- function(seurat_object, feature_num = 200,clustering_method = "ward.D", clustering_dist = "euclidean" ){
  ## process data
  DefaultAssay(seurat_object) <- "RNA"
  seurat_object <- FindVariableFeatures(seurat_object,nfeatures = 3000)
  heatdata <- GetAssayData(seurat_object, slot = "scale.data", assay = "RNA")[seurat_object@assays$SCT@var.features[1:feature_num],]
  ## Cell Annotation
  cellcolor_1 <- c("#0091B2","#3B4992FF", "#EEA236FF" )
  cellcolor_2 <- c("#3B4952FF","#EEA236FF", "#BB0021FF")
  seurat_object$Batch2 <- factor(seurat_object$Batch2, levels = c("Treated NPC43 (24hr)","Treated NPC43 (48hr)","Untreated NPC43"))
  seurat_object$nonresponsive <- factor(seurat_object$nonresponsive, levels = c("Induced state","Native state","Nonresponsive state"))
  group_col_1 <- structure(cellcolor_1, names = as.character(levels(seurat_object$Batch2)))
  group_col_2 <- structure(cellcolor_2, names = as.character(levels(seurat_object$nonresponsive)))
  ann_df <- data.frame(row.names = colnames(seurat_object), 'Cell_behavior' = seurat_object$nonresponsive, 'Treatment' = seurat_object$Batch2)
  ann <- HeatmapAnnotation(df = ann_df,
                           col = list('Treatment' = group_col_1, 'Cell_behavior' = group_col_2),
                           show_annotation_name = FALSE,
                           annotation_name_gp = gpar(fontsize = 8),
                           simple_anno_size = unit(0.5, "cm")
  )
  heatdata = as.matrix(heatdata)
  #color_heat <- structure(pals::ocean.balance(length(0:max(heatdata))),names = 0:max(heatdata))
  ht <- Heatmap(
    heatdata,
    use_raster = TRUE,
    raster_quality = 1,
    show_row_dend = F,
    show_column_dend = F,
    column_title = NULL,
    column_title_gp = gpar(fontsize = 12),
    column_title_side = "bottom",
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    clustering_distance_rows = clustering_dist,
    clustering_method_rows = clustering_method,
    cluster_columns = TRUE,
    cluster_column_slices = FALSE,
    clustering_distance_columns = clustering_dist,
    clustering_method_columns = clustering_method,
    border = "grey",
    row_title = NULL,
    row_title_gp = gpar(fontsize = 12),
    row_title_side = "right",
    show_column_names = FALSE,
    show_row_names = FALSE,
    #col = circlize::colorRamp2(c(0, 3, 6), c("#377EB8", "white", "#E41A1C")),
    #col = color_heat,
    top_annotation  = ann,
    show_heatmap_legend = T,
    #height = unit(2, "cm")
  )
  ht
}
getTOP_tfs = function(rss, celltypes, n = 20){
  rssThisType <- sort(rss[, celltypes], decreasing = TRUE)
  thisRss <- data.frame(regulon = names(rssThisType), rank = seq_along(rssThisType), rss = rssThisType)
  return(thisRss$regulon[1:n])
}
