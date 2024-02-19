# Script for analysing CosMx Data using TileDB/Seurat
# Written by James Boot 19/12/2023

# Set working directory
setwd('/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/NanostringCosMX/Projects/Xenakis_Dorian/GC-TX-10712/')

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
library(factoextra)
library(skmeans)
library(Matrix)
library(cluster)

# Path to tiledb directory
tiledb_Dir <- 'Data/GCTX10712/805c0136-6cdd-4604-ad65-cbe15c843c08_TileDB'

# Read in SOMACollection
tiledb_scdataset <- tiledbsc::SOMACollection$new(uri = tiledb_Dir, 
                                                 verbose = TRUE)

# Check loading 
names(tiledb_scdataset$somas)
names(tiledb_scdataset$somas$RNA$members)

# Get raw counts 
counts <- tiledb_scdataset$somas$RNA$X$members$counts$to_matrix(batch_mode = TRUE) 
dim(counts)

# Get norm counts
norm <- tiledb_scdataset$somas$RNA_normalized_efe5a2bb.b9ea.461a.93ea.af21c287d5c8_1$X$members$data$to_matrix(batch_mode = TRUE)
dim(norm)

# Get cell meta data
metadata <- tiledb_scdataset$somas$RNA$obs$to_dataframe()
dim(metadata)
head(metadata)

# QC metrics are stored in the meta data - explore

# Change working directory to save plots
setwd('/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/NanostringCosMX/Projects/Xenakis_Dorian/GC-TX-10712/Analysis/analysis1/')

# How many FOVs?
length(unique(metadata$fov))

# nCounts_RNA in all data
plot <- ggplot(metadata, aes(x=nCount_RNA)) + 
  geom_histogram(binwidth = 10) +
  ggtitle('nCount_RNA in all cells, binwidth 10') +
  scale_y_continuous(labels = label_number()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        title = element_text(size = 14))
plot
ggsave(
  plot = plot,
  filename = paste0('nCount_RNA.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
) 

# nFeature_RNA
plot <- ggplot(metadata, aes(x=nFeature_RNA)) + 
  geom_histogram(binwidth = 10) +
  ggtitle('nFeature_RNA in all cells, binwidth 10') +
  scale_y_continuous(labels = label_number()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        title = element_text(size = 14))
plot
ggsave(
  plot = plot,
  filename = paste0('nFeature_RNA.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)  

# Dot plot of nFeature vs nCount RNA
plot <- ggplot(metadata, aes(x=nCount_RNA, y=nFeature_RNA)) +
  geom_point(size=1, alpha=0.2) +
  ggtitle('nCount_RNA vs nFeature_RNA in all cells') +
  scale_y_continuous(labels = label_number()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        title = element_text(size = 14))
plot
ggsave(
  plot = plot,
  filename = paste0('nCount_v_nFeature.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Dot plot of nFeature vs nCount RNA with QC flags
plot <- ggplot(metadata, aes(x=nCount_RNA, y=nFeature_RNA, color = qcFlagsCellComplex)) +
  geom_point(size=1, alpha=0.2) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  ggtitle('nCount_RNA vs nFeature_RNA in all cells') +
  scale_y_continuous(labels = label_number()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        title = element_text(size = 14))
plot 
ggsave(
  plot = plot,
  filename = paste0('nCount_v_nFeature_qc1.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Dot plot of nFeature vs nCount RNA with QC flags
plot <- ggplot(metadata, aes(x=nCount_RNA, y=nFeature_RNA, color = qcFlagsCellCounts)) +
  geom_point(size=1, alpha=0.2) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  ggtitle('nCount_RNA vs nFeature_RNA in all cells') +
  scale_y_continuous(labels = label_number()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        title = element_text(size = 14))
plot
ggsave(
  plot = plot,
  filename = paste0('nCount_v_nFeature_qc2.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Density plot of nFeature vs nCount RNA
plot <- ggplot(metadata, aes(x=nCount_RNA, y=nFeature_RNA)) +
  geom_bin2d(bins = 100) +
  ggtitle('nCount_RNA vs nFeature_RNA in all cells') +
  scale_y_continuous(labels = label_number()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        title = element_text(size = 14)) +
  scale_fill_binned(type = "viridis")
plot
ggsave(
  plot = plot,
  filename = paste0('nCount_v_nFeature_density.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Summarise QC metrics
# Cell Counts = % of cells above the minimum counts per cell threshold set in module parameters.
# Prop Neg = % of cells with negative probe count proportion less than the threshold set in module parameters.
# Cell Complex = % of cells with ratio of total counts to the number of detected genes greater than the value set in module parameters.
# Area = % of cells with cell area which is not designated as an outlier according to the threshold set in module parameters.
# FOV = % of FOVs that were not flagged in FOV QC, according to module parameters set by user. Mean counts / cell = 100
qcMetrics <- data.frame(
  Metric = c(
    'qcFlagsRNACounts',
    'qcFlagsCellCounts',
    'qcFlagsCellPropNeg',
    'qcFlagsCellComplex',
    'qcFlagsCellArea',
    'qcFlagsFOV'
  ),
  Percentage_Pass = c(
    (sum(metadata$qcFlagsRNACounts == 'Pass') / nrow(metadata)) * 100,
    (sum(metadata$qcFlagsCellCounts == 'Pass') / nrow(metadata)) * 100,
    (sum(metadata$qcFlagsCellPropNeg == 'Pass') / nrow(metadata)) * 100,
    (sum(metadata$qcFlagsCellComplex == 'Pass') / nrow(metadata)) * 100,
    (sum(metadata$qcFlagsCellArea == 'Pass') / nrow(metadata)) * 100,
    (sum(metadata$qcFlagsFOV == 'Pass') / nrow(metadata)) * 100
  )
)
write.csv(qcMetrics, 'qcMetrics.csv')

# Summarise for FOV level
qcFOV <- metadata %>%
  group_by(fov) %>%
  summarise(mean_counts_per_cell = mean(nCount_RNA),
            mean_features_per_cell = mean(nFeature_RNA))
write.csv(qcFOV, 'FOVqcMetrics.csv')

# Histogram of average nCounts per cell in FOVs 
plot <- ggplot(qcFOV, aes(x = mean_counts_per_cell)) + 
  geom_histogram(binwidth = 10) +
  ggtitle('mean_counts_per_cell in all FOVs, binwidth 10') +
  scale_y_continuous(labels = label_number()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        title = element_text(size = 14))
plot
ggsave(
  plot = plot,
  filename = paste0('mean_counts_per_cell_FOV.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Summarise QC metrics per cluster
qcCluster <- metadata %>%
  group_by(nn_46182286.0cc0.4c7d.8a07.5186d16a8422_1_cluster_cluster_6c3ffc3b.47da.4b10.a800.dcd0e04a167c_1) %>%
  count(qcFlagsCellCounts) %>%
  arrange(nn_46182286.0cc0.4c7d.8a07.5186d16a8422_1_cluster_cluster_6c3ffc3b.47da.4b10.a800.dcd0e04a167c_1)

# Total number of clusters
length(unique(qcCluster$nn_46182286.0cc0.4c7d.8a07.5186d16a8422_1_cluster_cluster_6c3ffc3b.47da.4b10.a800.dcd0e04a167c_1))
qcClusterFilt <- qcCluster %>%
  filter(qcFlagsCellCounts == 'Pass')

# Total number of clusters after removing cells
length(unique(qcClusterFilt$nn_46182286.0cc0.4c7d.8a07.5186d16a8422_1_cluster_cluster_6c3ffc3b.47da.4b10.a800.dcd0e04a167c_1))

# Extract UMAP 
umap_TileDB <- as.data.frame(tiledb_scdataset$somas$RNA$obsm$members$dimreduction_approximateumap$to_matrix())
test <- tiledb_scdataset$somas$RNA$obs$to_dataframe()
head(test)
umap_TileDB$nFeature_RNA <- test$nFeature_RNA
umap_TileDB$qcFlagsCellCounts <- test$qcFlagsCellCounts
umap_TileDB$qcFlagsCellComplex <- test$qcFlagsCellComplex
umap_TileDB$cluster <- factor(test$nn_46182286.0cc0.4c7d.8a07.5186d16a8422_1_cluster_cluster_6c3ffc3b.47da.4b10.a800.dcd0e04a167c_1,
                              levels = c(1:200))
colnames(umap_TileDB)[1:2] <- c('UMAP1', 'UMAP2')
head(umap_TileDB)

# Plot UMAP with nFeature
plot <-
  ggplot(
    umap_TileDB,
    aes(x = UMAP1,
        y = UMAP2,
        color = nFeature_RNA)
  ) +
  geom_point(size = 0.5, alpha = 0.5) +
  labs(x = 'UMAP1', y = 'UMAP2') +
  scale_color_gradient(low = "gray90", high = "blue") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )
plot
ggsave(
  plot = plot,
  filename = paste0('UMAP_nFeature.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)  

# Plot UMAP with qcFlagsCellCounts  
plot <-
  ggplot(
    umap_TileDB,
    aes(x = UMAP1,
        y = UMAP2,
        color = qcFlagsCellCounts)
  ) +
  geom_point(size = 0.5, alpha = 0.2) +
  labs(x = 'UMAP1', y = 'UMAP2') +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
plot
ggsave(
  plot = plot,
  filename = paste0('UMAP_qcFlagsCellCounts.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Plot UMAP with qcFlagsCellComplex
plot <-
  ggplot(
    umap_TileDB,
    aes(x = UMAP1,
        y = UMAP2,
        color = qcFlagsCellComplex)
  ) +
  geom_point(size = 0.5, alpha = 0.2) +
  labs(x = 'UMAP1', y = 'UMAP2') +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
plot
ggsave(
  plot = plot,
  filename = paste0('UMAP_qcFlagsCellComplex.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Plot UMAP with cluster
plot <-
  ggplot(
    umap_TileDB,
    aes(x = UMAP1,
        y = UMAP2,
        color = cluster)
  ) +
  geom_point(size = 0.5, alpha = 0.2) +
  labs(x = 'UMAP1', y = 'UMAP2') +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )
plot
ggsave(
  plot = plot,
  filename = paste0('UMAP_cluster.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Plot the cell locations in 2D space as in the tissue
cellCoords <- tiledb_scdataset$somas$RNA$obs$to_dataframe(
  attrs = c("x_FOV_px", "y_FOV_px", "x_slide_mm", "y_slide_mm", 
            "slide_ID_numeric", "Run_Tissue_name", "fov",
            'nn_46182286.0cc0.4c7d.8a07.5186d16a8422_1_cluster_cluster_6c3ffc3b.47da.4b10.a800.dcd0e04a167c_1'))
head(cellCoords)
head(metadata)

# Bind the metadata and cellCoords dataframes based on row names
cellCoords2 <- merge(cellCoords, 
                     metadata[, c('nCount_RNA',
                                  'nFeature_RNA',
                                  'cell_id',
                                  'complexity',
                                  'qcFlagsCellComplex',
                                  'qcFlagsCellCounts',
                                  'qcFlagsFOV')], 
                     by = 'row.names')
head(cellCoords2)
colnames(cellCoords2)[9] <- 'cluster'

# Plot clusters
plot <-
  ggplot(
    cellCoords2,
    aes(x = x_slide_mm, y = y_slide_mm, color = cluster)
  ) +
  geom_point(alpha = 0.05, size = 0.01) +
  facet_wrap(~ Run_Tissue_name) +
  coord_equal() +
  labs(title = "Cell coordinates in XY space, coloured by cluster") + 
  theme(legend.position = 'none')
plot
ggsave(
  plot = plot,
  filename = paste0('CellCoords_cluster.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Plot cell counts qc flag
plot <-
  ggplot(
    cellCoords2,
    aes(x = x_slide_mm, y = y_slide_mm, color = qcFlagsCellCounts)
  ) +
  geom_point(alpha = 0.05, size = 0.01) +
  facet_wrap(~ Run_Tissue_name) +
  coord_equal() +
  labs(title = "Cell coordinates in XY space, coloured by qcFlagsCellCounts") + 
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(legend.position = 'right')
plot
ggsave(
  plot = plot,
  filename = paste0('CellCoords_qcFlagsCellCounts.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Plot cell complexity score
plot <-
  ggplot(
    cellCoords2,
    aes(x = x_slide_mm, y = y_slide_mm, colour = complexity)
  ) +
  geom_point(alpha = 0.05, size = 0.01) +
  facet_wrap(~ Run_Tissue_name) +
  coord_equal() +
  labs(title = "Cell coordinates in XY space, coloured by complexity score") + 
  theme(legend.position = 'right') +
  scale_colour_gradient(low = 'grey', high = 'red')
plot
ggsave(
  plot = plot,
  filename = paste0('CellCoords_complexity.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Plot cell nFeatures
plot <-
  ggplot(
    cellCoords2,
    aes(x = x_slide_mm, y = y_slide_mm, colour = nFeature_RNA)
  ) +
  geom_point(alpha = 0.05, size = 0.01) +
  facet_wrap(~ Run_Tissue_name) +
  coord_equal() +
  labs(title = "Cell coordinates in XY space, coloured by nFeature_RNA") + 
  theme(legend.position = 'right') +
  scale_colour_gradient(low = 'grey', high = 'red')
plot
ggsave(
  plot = plot,
  filename = paste0('CellCoords_nFeature.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Extract transcript coords
transcriptCoords <- tiledb::tiledb_array(
  tiledb_scdataset$somas$RNA$obsm$members$transcriptCoords$uri,
  return_as="data.frame")[]
head(transcriptCoords)

# Lets look at an FOV of interest
slide <- 2
fov <- 20

slideName <-
  unique(cellCoords$Run_Tissue_name[cellCoords$slide_ID_numeric ==
                                      slide])

fovCoords <- cellCoords[cellCoords$slide_ID_numeric == slide &
                          cellCoords$fov == fov, ]
fovTranscriptCoords <-
  transcriptCoords[transcriptCoords$slideID == slide &
                     transcriptCoords$fov == fov, ]

targetCounts <- table(fovTranscriptCoords$target)
targets <- names(targetCounts[which(targetCounts >= 1000)])

fovTranscriptCoords <-
  fovTranscriptCoords[fovTranscriptCoords$target %in%
                        targets,]

plot <- ggplot(fovCoords, aes(x = x_FOV_px, y = y_FOV_px)) +
  geom_point(alpha = 0.6,
             size = 0.5,
             color = "black") +
  geom_point(
    data = fovTranscriptCoords,
    mapping = aes(x = x_FOV_px,
                  y = y_FOV_px,
                  color = target),
    size = 0.5,
    alpha = 0.6
  ) +
  theme_bw() +
  coord_equal() +
  guides(colour = guide_legend(override.aes = list(size = 1,
                                                   alpha = 1))) +
  labs(color = "RNA Target",
       title = paste0("RNA Transcripts in\n",
                      slideName, "\nFOV", fov))
ggsave(
  plot = plot,
  filename = paste0('FOV_Transcript_Example.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Here we try filtering the data and the recluster

# Identify cell names to keep
passingCells <- row.names(metadata)[metadata$qcFlagsCellCounts == 'Pass']

# Subset the normalised counts matrix
normFilt <- norm[, passingCells]

# Create seurat object from normalised counts
seuratObj <- CreateSeuratObject(counts = normFilt, assay = 'CosMx')

# Find variable features and scale
seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 1000)
seuratObj[["CosMx"]]$data <- log(seuratObj[["CosMx"]]$counts + 1)
seuratObj <- ScaleData(seuratObj, block.size = 50000)

# Run PCA
seuratObj <- RunPCA(seuratObj, features = VariableFeatures(object = seuratObj))
ElbowPlot(seuratObj)

# Cluster
seuratObj <- FindNeighbors(seuratObj, dims = 1:5)
seuratObj <- FindClusters(seuratObj, resolution = 0.5)

# UMAP
seuratObj <- RunUMAP(seuratObj, dims = 1:5)
DimPlot(seuratObj, reduction = "umap")
ggsave(plot = last_plot(),
       filename = 'UMAP_Filt_clusters.tiff',
       height = 8,
       width = 8,
       dpi = 300,
       units = 'in')

# Extract meta data from seurat object
seuratMeta <- seuratObj@meta.data

# Add new clusters to coords
# Bind the metadata and cellCoords dataframes based on row names
row.names(cellCoords2) <- cellCoords2$Row.names
cellCoords2$Row.names <- NULL
cellCoords3 <- merge(cellCoords2, 
                     seuratMeta[, c('CosMx_snn_res.0.5', 'seurat_clusters')], 
                     by = 'row.names')
head(cellCoords3)

# Plot new clusters on tissue
plot <-
  ggplot(
    cellCoords3,
    aes(x = x_slide_mm, y = y_slide_mm, color = seurat_clusters)
  ) +
  geom_point(alpha = 0.05, size = 0.01) +
  facet_wrap(~ Run_Tissue_name) +
  coord_equal() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "Cell coordinates in XY space, coloured by seurat_clusters") + 
  theme(legend.position = 'right')
plot
ggsave(
  plot = plot,
  filename = paste0('CellCoords_NEWcluster.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Find markers
clust.markers <- FindAllMarkers(seuratObj, only.pos = TRUE)
clust.markers.filt <- clust.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) 
write.csv(clust.markers.filt, 'newClusters_markers.csv')  

clust.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
write.csv(top10, 'newClusters_markersTop10.csv')  




