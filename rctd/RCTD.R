# Script for Seurat RCTD deconvolution of spatial Visium data
# Uses outputs from CARD.R

# Packages
library(tidyverse)
library(Seurat)
library(spacexr)

# Load data setup seurat obj
scDat <- readRDS("scRef/Ctrl_scDat.RDS")
scMeta <- readRDS('scRef/Ctrl_scMeta.RDS')
sc_obj <- CreateSeuratObject(scDat)

# Remove cell types with less than 25 cells
typeRemove <- scMeta %>%
  count(CellType) %>%
  filter(n < 25) %>%
  pull(CellType)

barcodesRemove <- scMeta$Detected[scMeta$CellType %in% typeRemove]

scMeta <- scMeta[!scMeta$Detected %in% barcodesRemove,]
scDat <- scDat[!colnames(scDat) %in% barcodesRemove,]

# Extract information to pass to the RCTD Reference function
celltype <- as.factor(scMeta$CellType)
names(celltype) <- scMeta$Detected
nUMI <- sc_obj$nCount_RNA
names(nUMI) <- colnames(sc_obj)
reference <- Reference(scDat, celltype, nUMI)

# Load spatial data
spatialDatList <- readRDS('spatial/spatialDatList.RDS')
spatialLocsList <- readRDS('spatial/spatialLocsList.RDS')

# Create folder for outputs
outs <- 'RCTD_outs'
if (!dir.exists(outs)) {
  dir.create(outs)
}

# Loop through all samples for analysis
for (sample in names(spatialDatList)) {
  
  # Set up query with the RCTD function Spatial RNA
  spCounts <- spatialDatList[[sample]]
  coords <- spatialLocsList[[sample]][, c('x', 'y')]
  query <- SpatialRNA(coords, spCounts, colSums(spCounts))
  
  # Create RCTD
  RCTD <- create.RCTD(query, reference, max_cores = 2)
  
  # Run RTCD
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  
  # Save object
  saveRDS(RCTD, file = paste0(outs, '/', sample, '_RCTD.RDS'))
  
  # Save results
  write.csv(RCTD@results$results_df,
            file = paste0(outs, '/', sample, '_RCTD.csv'))
  
}





