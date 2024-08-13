# Script 2 - pre-processing of CosMx data
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

# Load meta data and normalised counts
metadata <- readRDS(file = paste0(analysis.dir, '/cellMeta.RDS'))
norm <- readRDS(file = paste0(analysis.dir, '/normCounts.RDS'))

# Filter data based on qcFlags

# Identify cell names to keep
passingCells <-
  row.names(metadata)[metadata$qcFlagsCellCounts == 'Pass' &
                        metadata$qcFlagsRNACounts == 'Pass' &
                        metadata$qcFlagsCellPropNeg == 'Pass' &
                        metadata$qcFlagsCellComplex == 'Pass' &
                        metadata$qcFlagsCellArea == 'Pass']

# Subset the normalised counts matrix
normFilt <- norm[, passingCells]

# Subset the metadata 
metadataFilt <- metadata[passingCells, ]

# Save normFilt as RDS object for later use
saveRDS(normFilt, file = paste0(analysis.dir, '/normCountsQCfilt.RDS'))

# Create seurat object from normalised counts
seuratObj <- CreateSeuratObject(counts = normFilt, assay = 'CosMx')

# Find variable features and scale
seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 1000)
seuratObj[["CosMx"]]$data <- log(seuratObj[["CosMx"]]$counts + 1)
seuratObj <- ScaleData(seuratObj, block.size = 50000)

# Run PCA
seuratObj <- RunPCA(seuratObj, features = VariableFeatures(object = seuratObj))
ElbowPlot(seuratObj, ndims = 50)

# Cluster
seuratObj <- FindNeighbors(seuratObj, dims = 1:20)
seuratObj <- FindClusters(seuratObj, resolution = 1)

# UMAP
seuratObj <- RunUMAP(seuratObj, dims = 1:20)
plt <- DimPlot(seuratObj, reduction = "umap")
ggsave(plot = plt,
       filename = paste0(analysis.dir, '/UMAP2_clusters.tiff'),
       height = 8,
       width = 8,
       dpi = 300,
       units = 'in')

# Add new clusters to metadataFilt
metadataFilt <- merge(
  metadataFilt,
  data.frame(
    row.names = rownames(seuratObj@meta.data),
    seurat_clusters = seuratObj@meta.data[, c('seurat_clusters')]
  ),
  by = 'row.names'
)
rownames(metadataFilt) <- metadataFilt$Row.names
metadataFilt$Row.names <- NULL

# Plot new clusters on tissue
plt <-
  ggplot(
    metadataFilt,
    aes(x = x_slide_mm, y = y_slide_mm, color = seurat_clusters)
  ) +
  geom_point(alpha = 0.05, size = 0.01) +
  facet_wrap(~ Run_Tissue_name) +
  coord_equal() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "Cell coordinates in XY space, coloured by seurat_clusters") + 
  theme(legend.position = 'right')

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/Spatial_clusters.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Identify hybrid cells
# Score cells for signatures
# Cells scoring highly for multiple signatures could be hybrids
# Sample contains, immune, tumour, endothelial,

# Immune sig
ImmuneSig <- list(
  c(
    'CD68',
    'CD163',
    'AIF1',
    'C1QC',
    'C1QA',
    'C1QB',
    'CD8A',
    'CD3G',
    'GZMA',
    'EOMES',
    'IL7',
    'IL7R',
    'TNFRSF9',
    'ICOS',
    'PTPRC'
  )
)

# Tumour sig
TumourSig <- list(
  c(
    'PIK3CA-intron-1',
    'PDGFRA-intron-1',
    'MYC-intron-2',
    'MET-intron-1',
    'MDM2-intron-9',
    'EGFR-intron-11',
    'CDK4-intron-1',
    'PDGFRA',
    'MYC',
    'MET',
    'EGFR',
    'GFAP'
  )
)

# Endothelial/Fibroblast sig
EndoFibroSig <- list(c(
  'CDH5',
  'CLEC14A',
  'ESAM',
  'PECAM1',
  'COL1A1',
  'COL3A1',
  'COL1A2',
  'COL4A2'
))

# Score all cells for all signatures 
seuratObj <-
  AddModuleScore(seuratObj,
                 features = ImmuneSig,
                 ctrl = 5,
                 name = 'ImmuneSig')

seuratObj <-
  AddModuleScore(seuratObj,
                 features = TumourSig,
                 ctrl = 5,
                 name = 'TumourSig')

seuratObj <-
  AddModuleScore(seuratObj,
                 features = EndoFibroSig,
                 ctrl = 5,
                 name = 'EndoFibroSig')

# Plot signature scores
plt <- VlnPlot(seuratObj, features = 'ImmuneSig1', pt.size = 0)
ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/Vln_cluster_immune.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

plt <- VlnPlot(seuratObj, features = 'TumourSig1', pt.size = 0)
ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/Vln_cluster_tumour.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

plt <- VlnPlot(seuratObj, features = 'EndoFibroSig1', pt.size = 0)
ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/Vln_cluster_EndoFibro.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Identify hybrids by finding cells with above mean signature scores in multiple signatures
sigScores <- seuratObj@meta.data[, c('ImmuneSig1',
                                     'TumourSig1',
                                     'EndoFibroSig1')]

sigScores <- sigScores %>%
  mutate(ImmuneSigStat = case_when(ImmuneSig1 > mean(ImmuneSig1) ~ 1,
                                   ImmuneSig1 < mean(ImmuneSig1) ~ 0)) %>%
  mutate(TumourSigStat = case_when(TumourSig1 > mean(TumourSig1) ~ 1,
                                   TumourSig1 < mean(TumourSig1) ~ 0)) %>%
  mutate(EndoFibroSig1 = case_when(EndoFibroSig1 > mean(EndoFibroSig1) ~ 1,
                                   EndoFibroSig1 < mean(EndoFibroSig1) ~ 0)) %>%
  mutate(HybridStat = rowSums(across(c(ImmuneSigStat, TumourSigStat, EndoFibroSig1)))) %>%
  mutate(HybridStatVerb = case_when(HybridStat > 1 ~ 'Hybrid',
                                    HybridStat <= 1 ~ 'Single'))

hybridSummary <- sigScores %>%
  group_by(HybridStat) %>%
  summarise(n = n())

# Save hybrids cell names to file
saveRDS(sigScores, file = paste0(analysis.dir, '/hybrids.RDS'))

# Add hybrid status to seurat object
seuratObj <-
  AddMetaData(
    seuratObj,
    data.frame(
      row.names = rownames(sigScores),
      HybridStatus = sigScores$HybridStatVerb
    )
  )

# Visualise hybrid status on UMAP
plt <- DimPlot(seuratObj, reduction = "umap", group.by = 'HybridStatus')
ggsave(plot = plt,
       filename = paste0(analysis.dir, '/UMAP2_HybridStatus.tiff'),
       height = 8,
       width = 8,
       dpi = 300,
       units = 'in')

# Plot hybrids on spatial plot
# Add hybrids to metadataFilt
metadataFilt <- merge(
  metadataFilt,
  data.frame(
    row.names = rownames(sigScores),
    HybridStatus = sigScores$HybridStatVerb
  ),
  by = 'row.names'
)
rownames(metadataFilt) <- metadataFilt$Row.names
metadataFilt$Row.names <- NULL

# Plot hybrids on tissue
plt <-
  ggplot(
    metadataFilt,
    aes(x = x_slide_mm, y = y_slide_mm, color = HybridStatus)
  ) +
  geom_point(alpha = 0.05, size = 0.01) +
  facet_wrap(~ Run_Tissue_name) +
  coord_equal() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "Cell coordinates in XY space, coloured by HybridStatus") + 
  theme(legend.position = 'right')

ggsave(
  plot = plt,
  filename = paste0(analysis.dir, '/Spatial_HybridStatus.tiff'),
  height = 10,
  width = 10,
  units = 'in',
  dpi = 300
)

# Save meta data
saveRDS(metadataFilt, file = paste0(analysis.dir, '/cellMetaQCfilt.RDS'))



