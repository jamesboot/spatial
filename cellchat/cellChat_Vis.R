# This is a script for creating specific CellChat visualisations for Milner et al.

# Load libraries
library(CellChat)
library(patchwork)
library(Seurat)

# Locate cellchat object RDS files 
inFiles <- list.files(
  path = 'cellChat_outs',
  pattern = '*.rds',
  recursive = T,
  full.names = T
)

# Create vector of sample names
samples <- gsub(
  '_cellchat.rds',
  '',
  basename(inFiles)
)

# Pathways to visualise
pathways <- c('Opioid', 'Tachykinins')

# Run through all samples in for loop
for (path in pathways) {
  for (x in 1:length(samples)) {
    obj <- readRDS(inFiles[x])
    if (path %in% obj@netP$pathways) {
      pdf(
        file = paste0(
          'cellChat_outs/',
          samples[x],
          '/',
          samples[x],
          '_',
          path,
          '_net.pdf'
        ),
        height = 5,
        width = 5
      )
      print(
        netVisual_aggregate(
          obj,
          signaling = path,
          layout = "spatial",
          edge.width.max = 2,
          vertex.size.max = 1,
          alpha.image = 0.2,
          vertex.label.cex = 3.5
        )
      )
      dev.off()
    } else {
      next
    }
  }
}

