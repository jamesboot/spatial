# Script for CARD deconvolution of spatial Visium data

# Packages
library(tidyverse)
library(Seurat)
library(CARD)

# Locate spatial data files
spatialFiles <- list.files(
  path = 'spatial',
  pattern = '*.h5',
  full.names = T,
  recursive = T
)

# Get sample names of h5 files
samples1 <- sub("_.*$",
                "",
                list.files(
                  path = 'spatial',
                  pattern = '*.h5',
                  recursive = T
                ))

# Read in spatial data
spatialDatList <- lapply(spatialFiles, function(x) {
  Read10X_h5(x)
})
names(spatialDatList) <- samples1

# Locate spatial locations files
spatialLocs <- list.files(
  path = 'spatial',
  pattern = 'tissue_positions.csv',
  full.names = T,
  recursive = T
)

# Get sample names of locs files
samples2 <-
  gsub('_spatial/spatial/tissue_positions.csv', '', spatialLocs)
samples2 <- gsub('spatial/', '', samples2)

# Load spatial locations
spatialLocsList <- lapply(spatialLocs, function(x) {
  read.csv(x)
})

names(spatialLocsList) <- samples2

# Edit spatial locs dataframes
for (x in names(spatialLocsList)) {
  rownames(spatialLocsList[[x]]) <- spatialLocsList[[x]]$barcode
  colnames(spatialLocsList[[x]]) <-
    c("barcode",
      "in_tissue",
      "x",
      "y",
      "pxl_row_in_fullres",
      "pxl_col_in_fullres")
}

# Subset spatial down to filtered spots
selectSpots <- readRDS('spatial/selectedSpots.RDS')
for (x in names(spatialDatList)) {
  tmp <- selectSpots[selectSpots$`Sample name` == x,]
  colsToKeep <-
    colnames(spatialDatList[[x]])[colnames(spatialDatList[[x]]) %in% tmp$`Cell name`]
  spatialDatList[[x]] <- spatialDatList[[x]][, colsToKeep]
  spatialLocsList[[x]] <- spatialLocsList[[x]][colsToKeep,]
  rm(tmp)
  rm(colsToKeep)
}

# Read in single cell dat
scDat <-
  as.matrix(read.delim(
    gzfile('scRef/GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt.gz')
  ))
scMeta <-
  read.delim(gzfile('scRef/GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz'))

# Correct cellnames in meta to match counts mat
scMeta$Detected <- gsub(':', '.', scMeta$Detected)
scMeta$Detected <- gsub('/', '.', scMeta$Detected)

# Check that cell names match
if (sum(colnames(scDat) %in% scMeta$Detected) == nrow(scMeta)) {
  message('SUCCESS: All cell names match.')
} else if (sum(colnames(scDat) %in% scMeta$Detected) != nrow(scMeta)) {
  message('ERROR: All cell names do not match.')
}

# Only merge cell types e.g. neuron1, neuron2 = neuron
# Define cell type groups
Neurons <- c("Neuron1", "Neuron2", "Neuron3", "Neuron4", "Neuron5")
Oligo <-
  c("Oligo2", "Oligo6", "Oligo1", "Oligo4", "Oligo3", "Oligo5", "COPs", "ImOlGs", "OPCs")
Astros <- c("Astrocytes", "Astrocytes2")
Endo <- c("Endothelial_cells1", "Endothelial_cells2", "Pericytes", "Vasc_smooth_muscle")
Micro_Macros <- c("Microglia_Macrophages", "Macrophages")

scMeta <- scMeta %>%
  mutate(
    CellType = case_when(
      Celltypes %in% Neurons ~ "Neuron",
      Celltypes %in% Oligo ~ "Oligodendrocyte",
      Celltypes %in% Astros ~ "Astrocyte",
      Celltypes %in% Endo ~ "Endothelial_cell",
      Celltypes %in% Micro_Macros ~ "Microglia_Macrophage",
      Celltypes %in% c('Immune_cells') ~ "Immune_cells"
    )
  )

rownames(scMeta) <- scMeta$Detected

# Now subset down the matrix to selected cells from meta
scDat <- scDat[, scMeta$Detected]

# Save all objects for future use
saveRDS(scDat, file = 'scRef/Ctrl_scDat.RDS')
saveRDS(scMeta, file = 'scRef/Ctrl_scMeta.RDS')
saveRDS(spatialDatList, file = 'spatial/spatialDatList.RDS')
saveRDS(spatialLocsList, file = 'spatial/spatialLocsList.RDS')

# Create folder for outputs
outs <- 'CARD_outs'
if (!dir.exists(outs)) {
  dir.create(outs)
}

# Run through all samples
for (sample in names(spatialDatList)) {
  # Create CARD object
  CARD_obj = createCARDObject(
    sc_count = scDat,
    sc_meta = scMeta,
    spatial_count = spatialDatList[[sample]],
    spatial_location = spatialLocsList[[sample]],
    ct.varname = "CellType",
    ct.select = unique(scMeta$CellType),
    sample.varname = NULL,
    minCountGene = 100,
    minCountSpot = 5
  )
  
  # Deconvolve
  CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)
  
  # Visualise
  colors <- c(
    '#e6194b',
    '#3cb44b',
    '#ffe119',
    '#4363d8',
    '#f58231',
    '#911eb4',
    '#42d4f4',
    '#f032e6',
    '#bcf60c',
    '#fabebe',
    '#008080'
  )
  
  # Export proportions
  write.csv(CARD_obj@Proportion_CARD,
            file = paste0(outs, '/', sample, '.csv'))
  
  # Visualise
  p1 <- CARD.visualize.pie(
    proportion = CARD_obj@Proportion_CARD,
    spatial_location = CARD_obj@spatial_location,
    radius = 0.52,
    colors = colors
  )
  pdf(
    file = paste0(outs, '/', sample, '_CARD.pdf'),
    height = 10,
    width = 10
  )
  print(p1)
  dev.off()
}
