# Script 3 - post-processing of CosMx data
# Written by James Boot 17/05/2024

# Set working directory
setwd('/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/NanostringCosMX/Projects/Noorani_Imran/GC-IN-10692/Analysis/')
analysis.dir <- '/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/NanostringCosMX/Projects/Noorani_Imran/GC-IN-10692/Analysis/analysis6'

# Create folder for results to go in
if (!dir.exists(analysis.dir)) {
  dir.create(analysis.dir)
} 

# Load packages
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(matrixStats)
library(patchwork)
library(pheatmap)
library(Seurat)
library(RColorBrewer)
library(reshape2)
library(tiledb)
library(tiledbsc)
library(ggbreak)
library(dplyr)
library(scales)
library(ComplexHeatmap)
library(openxlsx)
library(EnhancedVolcano)

# Load filtered meta data and normalised counts
metadata <- readRDS(file = paste0(analysis.dir, '/cellMetaQCfilt.RDS'))
norm <- readRDS(file = paste0(analysis.dir, '/normCountsQCfilt.RDS'))
hybrids <- readRDS(file = paste0(analysis.dir, '/hybrids.RDS'))

# Create Seurat object
seuratObj <- CreateSeuratObject(counts = norm, assay = 'CosMx')

# Add hybrid status to seurat object
seuratObj <-
  AddMetaData(
    seuratObj,
    data.frame(
      row.names = rownames(hybrids),
      HybridStatus = hybrids$HybridStatVerb
    )
  )

# Get number of hybrid cells and total number of cells
sum(seuratObj$HybridStatus == 'Hybrid')
length(seuratObj$HybridStatus)

# Remove hybrids 
seuratObj <- subset(x = seuratObj, subset = HybridStatus == "Single")

# Find variable features and scale
seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 1000)
seuratObj[["CosMx"]]$data <- log(seuratObj[["CosMx"]]$counts + 1)
seuratObj <- ScaleData(seuratObj, block.size = 50000)

# Run PCA
seuratObj <- RunPCA(seuratObj, features = VariableFeatures(object = seuratObj))
ElbowPlot(seuratObj)

# Cluster
seuratObj <- FindNeighbors(seuratObj, dims = 1:20)
seuratObj <- FindClusters(seuratObj, resolution = 1)

# UMAP
seuratObj <- RunUMAP(seuratObj, dims = 1:20)
plt <- DimPlot(seuratObj, reduction = "umap")
ggsave(plot = plt,
       filename = paste0(analysis.dir, '/UMAP3_clusters.tiff'),
       height = 8,
       width = 8,
       dpi = 300,
       units = 'in')

# Add new clusters to metadata
metadata <- merge(
  metadata,
  data.frame(
    row.names = rownames(seuratObj@meta.data),
    seurat_clusters_final = seuratObj@meta.data[, c('seurat_clusters')]
  ),
  by = 'row.names'
)
rownames(metadata) <- metadata$Row.names
metadata$Row.names <- NULL

# Plot new clusters on tissue
plt <-
  ggplot(
    metadata,
    aes(x = x_slide_mm, y = y_slide_mm, color = seurat_clusters_final)
  ) +
  geom_point(alpha = 0.05, size = 0.01) +
  facet_wrap(~ Run_Tissue_name) +
  coord_equal() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "Cell coordinates in XY space, coloured by seurat_clusters_final") + 
  theme(legend.position = 'right')

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/Spatial_clusters_noHybrid.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Find markers of clusters 
clust.markers <- FindAllMarkers(seuratObj, only.pos = TRUE)
clust.markers.filt <- clust.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1 &
                  p_val_adj < 0.01)   %>%
  arrange(desc(avg_log2FC), .by_group = T)

write.csv(clust.markers.filt, file = paste0(analysis.dir, '/final_cluster_markers.csv'))

# Annotate cells using custom scType script

# scType
library(HGNChelper)

# Load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# Load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file - use custom file modified specifically for this project
db_ = "Copy of ScTypeDB_full.xlsx";
tissue = "Custom" 

# Prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# Get cell-type by cell matrix
es.max = sctype_score(
  scRNAseqData = seuratObj[['CosMx']]$scale.data,
  scaled = TRUE,
  gs = gs_list$gs_positive,
  gs2 = gs_list$gs_negative
) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data 

# Merge by cluster
cL_resutls = do.call("rbind", lapply(unique(seuratObj@meta.data$seurat_clusters), function(cl) {
  es.max.cl = sort(rowSums(es.max[, rownames(seuratObj@meta.data[seuratObj@meta.data$seurat_clusters ==
                                                                   cl, ])]), decreasing = !0)
  head(data.frame(
    cluster = cl,
    type = names(es.max.cl),
    scores = es.max.cl,
    ncells = sum(seuratObj@meta.data$seurat_clusters == cl)
  ),
  10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells /
                     4] = "Unknown"
print(sctype_scores[, 1:3])

# Plot results - results added as customclassif
seuratObj@meta.data$customclassif = ""
for (j in unique(sctype_scores$cluster)) {
  cl_type = sctype_scores[sctype_scores$cluster == j, ]
  
  seuratObj@meta.data$customclassif[seuratObj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

# Plot UMAP with classifications
plt <-
  DimPlot(
    seuratObj,
    reduction = "umap",
    label = TRUE,
    repel = TRUE,
    group.by = 'customclassif'
  )

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/UMAP3_scType1.tiff'),
  height = 8,
  width = 12,
  dpi = 300,
  units = 'in'
)

plt <-
  DimPlot(
    seuratObj,
    reduction = "umap",
    label = TRUE,
    repel = TRUE,
    split.by = 'customclassif'
  ) +
  theme(strip.text.x = element_text(size = 6, face = "bold"))

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/UMAP3_scType2.tiff'),
  height = 4,
  width = 16,
  dpi = 300,
  units = 'in'
)

saveRDS(seuratObj, file = paste0(analysis.dir, '/scTypeSeuratObj.RDS'))

# Manual annotation based on marker expression & scType
manualClass <- c()
manualClass[seuratObj$seurat_clusters %in% c(0)] <- 'Glia/Schwann'
manualClass[seuratObj$seurat_clusters %in% c(1,2,3,6,7)] <- 'Tumour cells'
manualClass[seuratObj$seurat_clusters %in% c(4,5,18)] <- 'Macrophages'
manualClass[seuratObj$seurat_clusters %in% c(8)] <- 'Immature neurons'
manualClass[seuratObj$seurat_clusters %in% c(9,10,19)] <- 'Glia'
manualClass[seuratObj$seurat_clusters %in% c(11)] <- 'Microglia'
manualClass[seuratObj$seurat_clusters %in% c(12,16)] <- 'T cells'
manualClass[seuratObj$seurat_clusters %in% c(13)] <- 'Mesenchymal/glia'
manualClass[seuratObj$seurat_clusters %in% c(14)] <- 'Neurons/glia'
manualClass[seuratObj$seurat_clusters %in% c(15)] <- 'Neurons'
manualClass[seuratObj$seurat_clusters %in% c(17)] <- 'Astrocytes'
unique(manualClass)
seuratObj$Manual_Anno <- manualClass

# UMAP
plt <-
  DimPlot(seuratObj, reduction = "umap", group.by = 'Manual_Anno')
ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/UMAP3_manualAnno1.tiff'),
  height = 8,
  width = 8,
  dpi = 300,
  units = 'in'
)

plt <-
  DimPlot(
    seuratObj,
    reduction = "umap",
    label = TRUE,
    repel = TRUE,
    split.by = 'Manual_Anno'
  ) +
  theme(strip.text.x = element_text(size = 6, face = "bold"))

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/UMAP3_manualAnno2.tiff'),
  height = 4,
  width = 16,
  dpi = 300,
  units = 'in'
)

# Add manual anno to meta 
metadata <- merge(
  metadata,
  data.frame(
    row.names = rownames(seuratObj@meta.data),
    Manual_Anno = seuratObj@meta.data[, c('Manual_Anno')]
  ),
  by = 'row.names'
)
rownames(metadata) <- metadata$Row.names
metadata$Row.names <- NULL

# Plot new clusters on tissue
plt <-
  ggplot(
    metadata,
    aes(x = x_slide_mm, y = y_slide_mm, color = Manual_Anno)
  ) +
  geom_point(alpha = 0.05, size = 0.01) +
  facet_wrap(~ Run_Tissue_name) +
  coord_equal() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "Cell coordinates in XY space, coloured by Manual_Anno") + 
  theme(legend.position = 'right')

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/Spatial_manualAnno.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Look at EGFR Expression distribution

# Subset down to tumour cells
tumourObj <- subset(x = seuratObj, subset = Manual_Anno == "Tumour cells")

# Extract EGFR expression
egfrExp <- data.frame(row.names = colnames( tumourObj[["CosMx"]]$counts),
                      EGFR = tumourObj[["CosMx"]]$counts['EGFR', ])

plt <- ggplot(egfrExp, aes(EGFR)) +
  geom_histogram()

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/EGFR_distrbtn.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Use threshold of 3
egfrExp <- egfrExp %>%
  mutate(EGFRStat = case_when(EGFR >= 3 ~ 'EGFR+',
                              EGFR < 3 ~ 'EGFR-'))

# Add EGFR status to meta data of tumour object
tumourObj <- AddMetaData(tumourObj,
                         data.frame(row.names = rownames(egfrExp),
                                    EGFR_Status = egfrExp$EGFRStat))

# Add FOV to meta data
tumourObj <- AddMetaData(tumourObj,
                         data.frame(row.names = rownames(metadata),
                                    FOV = paste0(metadata$slide_ID_numeric, '_', metadata$fov)))

# View EGFR status on umap
plt <- DimPlot(tumourObj, reduction = "umap", group.by = 'EGFR_Status')
ggsave(plot = plt,
       filename = paste0(analysis.dir, '/UMAP3_EGFR_Status.tiff'),
       height = 8,
       width = 8,
       dpi = 300,
       units = 'in')

# Compare EGFR+ Tumour cells with EGFR- using single cells
Idents(tumourObj) <- 'EGFR_Status'
egfr_degs <- FindMarkers(object = tumourObj, 
                         ident.1 = 'EGFR+',
                         ident.2 = 'EGFR-',
                         only.pos = F,
                         test.use = "wilcox")


# # Compare EGFR+ Tumour cells with EGFR- psuedobulk
# tumourBulk <-
#   AggregateExpression(
#     tumourObj,
#     assays = "CosMx",
#     return.seurat = T,
#     group.by = c("FOV", "EGFR_Status")
#   )
# 
# tail(Cells(tumourBulk))
# 
# Idents(tumourBulk) <- "EGFR_Status"
# egfr_degs <- FindMarkers(object = tumourBulk, 
#                             ident.1 = "EGFR+", 
#                             ident.2 = "EGFR-",
#                             test.use = "wilcox")
# head(egfr_degs, n = 15)
# tail(egfr_degs, n = 15)

write.csv(egfr_degs, file = paste0(analysis.dir, '/EGFR_DEGs.csv'))

# Prepare custom colours
egfr_degs$colour <- NA
egfr_degs$desc <- NA
for (i in 1:nrow(egfr_degs)) {
  if (is.na(egfr_degs$p_val_adj[i])) {
    egfr_degs$colour[i] <- 'grey'
      egfr_degs$desc[i] <- 'NS'
  } else if (egfr_degs$p_val_adj[i] < 0.01 && egfr_degs$avg_log2FC[i] > 1) {
    egfr_degs$colour[i] <- 'red'
      egfr_degs$desc[i] <- 'Up-regulated'
  } else if (egfr_degs$p_val_adj[i] < 0.01 && egfr_degs$avg_log2FC[i] < -1) {
    egfr_degs$colour[i] <- 'blue'
      egfr_degs$desc[i] <- 'Down-regulated'
  } else {
    egfr_degs$colour[i] <- 'grey'
      egfr_degs$desc[i] <- 'NS'
  }
}

colours <- egfr_degs$colour
names(colours) <- egfr_degs$desc
egfr_degs$gene <- rownames(egfr_degs)

volc <- EnhancedVolcano(egfr_degs,
                        lab = egfr_degs$gene,
                        selectLab = egfr_degs$gene[1:20],
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        title = 'EGFR+ Tumour cells vs EGFR+ Tumour cells',
                        subtitle = bquote(italic('Pseudo-bulk')),
                        axisLabSize = 12,
                        pCutoff = 0.01,
                        FCcutoff = 1,
                        labSize = 3,
                        labCol = 'black',
                        labFace = 'bold',
                        boxedLabels = TRUE,
                        colAlpha = 4/5,
                        drawConnectors = TRUE,
                        colConnectors = 'black',
                        border = 'full',
                        borderColour = 'black',
                        borderWidth = 1,
                        legendLabSize = 10,
                        legendIconSize = 5,
                        colCustom = colours
)

ggsave(
  plot = volc,
  filename = paste0(analysis.dir, '/degs_EGFR.tiff'),
  height = 10,
  width = 10,
  dpi = 300,
  units = 'in'
)

# Investigate CHEK1 expression
seuratObj <- readRDS('analysis6/finalSeuratObj.RDS')

# Correlate CHEK1 and EGFR
corrDf <- data.frame(
  CHEK1 = seuratObj[['CosMx']]$counts['CHEK1',],
  EGFR = seuratObj[['CosMx']]$counts['EGFR',],
  EGFRcustom = seuratObj[['CosMx']]$counts['EGFR-intron-11',],
  PDGFRA = seuratObj[['CosMx']]$counts['PDGFRA',],
  PDGFRAcustom = seuratObj[['CosMx']]$counts['PDGFRA-intron-1',],
  MYC = seuratObj[['CosMx']]$counts['MYC',],
  MYCcustom = seuratObj[['CosMx']]$counts['MYC-intron-2',],
  MET = seuratObj[['CosMx']]$counts['MET',],
  METcustom = seuratObj[['CosMx']]$counts['MET-intron-1',],
  PIK3CAcustom = seuratObj[['CosMx']]$counts['PIK3CA-intron-1',],
  BGN = seuratObj[['CosMx']]$counts['BGN',],
  random1 = seuratObj[['CosMx']]$counts[sample(rownames(seuratObj[['CosMx']]$counts), 1),],
  random2 = seuratObj[['CosMx']]$counts[sample(rownames(seuratObj[['CosMx']]$counts), 1),],
  random3 = seuratObj[['CosMx']]$counts[sample(rownames(seuratObj[['CosMx']]$counts), 1),],
  random4 = seuratObj[['CosMx']]$counts[sample(rownames(seuratObj[['CosMx']]$counts), 1),],
  random5 = seuratObj[['CosMx']]$counts[sample(rownames(seuratObj[['CosMx']]$counts), 1),],
  random6 = seuratObj[['CosMx']]$counts[sample(rownames(seuratObj[['CosMx']]$counts), 1),],
  random7 = seuratObj[['CosMx']]$counts[sample(rownames(seuratObj[['CosMx']]$counts), 1),]
)

# For loop to go through all genes of interest
# Make dataframe for correlations
CorRes <- data.frame(Correlation = NA)
for (z in c(2:ncol(corrDf))) {
  gene <- colnames(corrDf)[z]
  # Plot correlation
  plt <- ggplot(corrDf, aes_string(x = gene, y = 'CHEK1')) +
    geom_point(alpha = 0.1)
  ggsave(
    plot = plt,
    filename = paste0(analysis.dir, '/corr_', gene, '.tiff'),
    height = 8,
    width = 8,
    dpi = 300,
    units = 'in'
  )
  # Calculate correlation
  CorRes[gene,  'Correlation'] <- cor(corrDf[, gene], corrDf[, 'CHEK1'])
}
write.csv(CorRes, paste0(analysis.dir,  '/CHEK1_correlations.csv'))

# End of analysis 

# Save seurat object
saveRDS(seuratObj, paste0(analysis.dir, '/finalSeuratObj.RDS'))

# Export data needed for InSituCor - run in InSituCor script
countMat <- seuratObj[["CosMx"]]$counts
saveRDS(countMat, paste0(analysis.dir, '/InSituCor_countMat.RDS'))
saveRDS(metadata, paste0(analysis.dir, '/InSituCor_meta.RDS'))
