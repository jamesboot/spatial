# Script 1 - QC of CosMx data
# Written by James Boot 17/05/2024

# Set working directory
setwd('/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/NanostringCosMX/Projects/Noorani_Imran/GC-IN-10692')
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

# Path to tiledb directory
tiledb_Dir <- 'Data/GCIN10692Exp122InitialtileDB/63c710fa-e431-4dcd-b852-87d780af305f_TileDB/'

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
norm <- tiledb_scdataset$somas$RNA_normalized_41009d84.5c3f.4c96.9003.c1e8ebfeccdd_1$X$members$data$to_matrix(batch_mode = TRUE)
dim(norm)

# Get cell meta data
metadata <- tiledb_scdataset$somas$RNA$obs$to_dataframe()
dim(metadata)
head(metadata)

# Save objects created
saveRDS(counts, file = paste0(analysis.dir, '/rawCounts.RDS'))
saveRDS(norm, file =  paste0(analysis.dir, '/normCounts.RDS'))
saveRDS(metadata, file =  paste0(analysis.dir, '/cellMeta.RDS'))

# QC metrics are stored in the meta data - explore

# nCounts_RNA in all data
plt <- ggplot(metadata, aes(x=nCount_RNA)) + 
  geom_histogram(binwidth = 10) +
  ggtitle('nCount_RNA in all cells, binwidth 10') +
  scale_y_continuous(labels = label_number()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        title = element_text(size = 14))
ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/QC_nCount.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
) 

# nFeature_RNA
plt <- ggplot(metadata, aes(x=nFeature_RNA)) + 
  geom_histogram(binwidth = 10) +
  ggtitle('nFeature_RNA in all cells, binwidth 10') +
  scale_y_continuous(labels = label_number()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        title = element_text(size = 14))
ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/QC_nFeature.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)  

# Dot plot of nFeature vs nCount RNA
plt <- ggplot(metadata, aes(x=nCount_RNA, y=nFeature_RNA)) +
  geom_point(size=1, alpha=0.2) +
  ggtitle('nCount_RNA vs nFeature_RNA in all cells') +
  scale_y_continuous(labels = label_number()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        title = element_text(size = 14))
ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/QC_nCount_v_nFeature1.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Dot plot of nFeature vs nCount RNA with QC flags
plt <- ggplot(metadata, aes(x=nCount_RNA, y=nFeature_RNA, color = qcFlagsCellComplex)) +
  geom_point(size=1, alpha=0.2) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  ggtitle('nCount_RNA vs nFeature_RNA in all cells') +
  scale_y_continuous(labels = label_number()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        title = element_text(size = 14))
ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/QC_nCount_v_nFeature2.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Dot plot of nFeature vs nCount RNA with QC flags
plt <- ggplot(metadata, aes(x=nCount_RNA, y=nFeature_RNA, color = qcFlagsCellCounts)) +
  geom_point(size=1, alpha=0.2) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  ggtitle('nCount_RNA vs nFeature_RNA in all cells') +
  scale_y_continuous(labels = label_number()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        title = element_text(size = 14))
ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/QC_nCount_v_nFeature3.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Density plot of nFeature vs nCount RNA
plt <- ggplot(metadata, aes(x=nCount_RNA, y=nFeature_RNA)) +
  geom_bin2d(bins = 100) +
  ggtitle('nCount_RNA vs nFeature_RNA in all cells') +
  scale_y_continuous(labels = label_number()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        title = element_text(size = 14)) +
  scale_fill_binned(type = "viridis")

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/QC_nCount_v_nFeature4.tiff'),
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
write.csv(qcMetrics, file = paste0(analysis.dir, '/QC_Metrics.csv'))

# Summarise for FOV level
qcFOV <- metadata %>%
  group_by(fov) %>%
  summarise(mean_counts_per_cell = mean(nCount_RNA),
            mean_features_per_cell = mean(nFeature_RNA))
write.csv(qcFOV, file = paste0(analysis.dir, '/QC_FOV_Metrics.csv'))

# Histogram of average nCounts per cell in FOVs 
plt <- ggplot(qcFOV, aes(x = mean_counts_per_cell)) + 
  geom_histogram(binwidth = 10) +
  ggtitle('mean_counts_per_cell in all FOVs, binwidth 10') +
  scale_y_continuous(labels = label_number()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        title = element_text(size = 14))

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/QC_nCounts_FOV.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
) 

# Extract UMAP 
umap_TileDB <- as.data.frame(tiledb_scdataset$somas$RNA$obsm$members$dimreduction_approximateumap$to_matrix())
umap_TileDB$nFeature_RNA <- metadata$nFeature_RNA
umap_TileDB$qcFlagsCellCounts <- metadata$qcFlagsCellCounts
umap_TileDB$qcFlagsCellComplex <- metadata$qcFlagsCellComplex
umap_TileDB$cluster <- factor(metadata$nn_71a26cad.df49.4509.87bf.d0ed68637f66_1_cluster_cluster_bde1f611.ef8b.40d0.b457.8a7801e1d016_1,
                              levels = c(1:200))
colnames(umap_TileDB)[1:2] <- c('UMAP1', 'UMAP2')
head(umap_TileDB)

# Plot UMAP with nFeature
plt <-
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

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/UMAP1_nFeature.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)  

# Plot UMAP with qcFlagsCellCounts  
plt <-
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

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/UMAP1_qc1.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Plot UMAP with qcFlagsCellComplex
plt <-
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

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/UMAP1_qc2.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Plot UMAP with cluster
plt <-
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

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/UMAP1_cluster.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Plot the cell locations in 2D space as in the tissue

# Plot clusters
colnames(metadata)[1] <- 'cluster'
plt <-
  ggplot(
    metadata,
    aes(x = x_slide_mm, y = y_slide_mm, color = cluster)
  ) +
  geom_point(alpha = 0.05, size = 0.01) +
  facet_wrap(~ Run_Tissue_name) +
  coord_equal() +
  labs(title = "Cell coordinates in XY space, coloured by cluster") + 
  theme(legend.position = 'none')

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/Spatial_cluster.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Plot cell counts qc flag
plt <-
  ggplot(
    metadata,
    aes(x = x_slide_mm, y = y_slide_mm, color = qcFlagsCellCounts)
  ) +
  geom_point(alpha = 0.05, size = 0.01) +
  facet_wrap(~ Run_Tissue_name) +
  coord_equal() +
  labs(title = "Cell coordinates in XY space, coloured by qcFlagsCellCounts") + 
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(legend.position = 'right')

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/Spatial_qc1.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Plot cell complexity score
plt <-
  ggplot(
    metadata,
    aes(x = x_slide_mm, y = y_slide_mm, colour = complexity)
  ) +
  geom_point(alpha = 0.05, size = 0.01) +
  facet_wrap(~ Run_Tissue_name) +
  coord_equal() +
  labs(title = "Cell coordinates in XY space, coloured by complexity score") + 
  theme(legend.position = 'right') +
  scale_colour_gradient(low = 'grey', high = 'red')

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/Spatial_qc2.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Plot cell nFeatures
plt <-
  ggplot(
    metadata,
    aes(x = x_slide_mm, y = y_slide_mm, colour = nFeature_RNA)
  ) +
  geom_point(alpha = 0.05, size = 0.01) +
  facet_wrap(~ Run_Tissue_name) +
  coord_equal() +
  labs(title = "Cell coordinates in XY space, coloured by nFeature_RNA") + 
  theme(legend.position = 'right') +
  scale_colour_gradient(low = 'grey', high = 'red')

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/Spatial_nFeature.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)






