# This is a script for running CellChat for Milner et al.

# Load libraries
library(CellChat)
library(patchwork)
library(Seurat)

options(stringsAsFactors = FALSE)
maxsize <- 4000*1024^2
options(future.globals.maxSize=maxsize)

# Set outdir
outdir <- 'cellChat_outs/'

# Locate all files
sampleDirs <- list.dirs(path = 'spatial', recursive = F)

# Create vector of sample names from sampleDirs
samples <- gsub('spatial/', '', sampleDirs)
samples <- gsub('_spatial', '', samples)

# Later we will need json files for each sample
jsonFiles <- list.files(
  path = 'spatial',
  pattern = '*.json',
  recursive = T,
  full.names = T
)

# For loop to process all samples and save cellchat RDS object for each
for (x in 11:length(samples)) {
  # Log
  message(paste('Starting sample:', samples[x]))
  
  # Set sample outdir
  sampledir <- file.path(outdir, samples[x])
  if (!dir.exists(sampledir)) {
    dir.create(sampledir)
  }
  
  # Load data
  obj <- Load10X_Spatial(
    data.dir = sampleDirs[x],
    filename = paste0(samples[x], '_filtered_feature_bc_matrix.h5'),
    assay = 'Spatial',
    filter.matrix = T,
    to.upper = F
  )
  
  # Add cluster meta data from Partek
  message(paste('Adding meta data for:', samples[x]))
  partekMeta <- read.delim('exported.txt')
  partekMetaFilt <- partekMeta[partekMeta$Sample.name == samples[x],]
  clustersMeta <- data.frame(
    row.names = partekMetaFilt$Cell.name,
    PartekCluster = as.factor(partekMetaFilt$Clusters..5comp.0.75res)
  )
  obj <- AddMetaData(obj, clustersMeta)
  
  # Filter to only cells
  obj <- subset(obj, cells = row.names(clustersMeta))
  
  # Pre-process visium data
  message(paste('Seurat processing sample:', samples[x]))
  obj <- SCTransform(obj, assay = "Spatial", verbose = T)
  obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
  obj <- FindClusters(obj, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
  
  # Prepare input data for CelChat analysis
  message(paste('Preparing cellchat input for sample:', samples[x]))
  data.input = Seurat::GetAssayData(obj, layer = "data", assay = "SCT")
  
  # Define the meta data:
  # A column named `slices` should be provided for spatial transcriptomics analysis,
  # `slices` is useful for analyzing cell-cell communication by aggregating multiple slices/replicates.
  # For comparison analysis across different conditions, users still need to create a CellChat object separately for each condition.
  # Ensure Idents of seurat object is set to Partek clusters
  Idents(object = obj) <- "PartekCluster"
  # Manually create a dataframe consisting of the cell labels
  meta = data.frame(
    labels = Seurat::Idents(obj),
    slices = "slice1",
    row.names = names(Seurat::Idents(obj))
  )
  meta$slices <- factor(meta$slices)
  unique(meta$labels)
  
  # Get tissue coordinates
  spatial.locs = Seurat::GetTissueCoordinates(obj,
                                              scale = NULL,
                                              cols = c("imagerow", "imagecol"))[c('x', 'y')]
  
  
  # Calculate scale factors
  scalefactors = jsonlite::fromJSON(txt = jsonFiles[x])
  # The theoretical spot size (um) in 10X Visium
  spot.size = 65
  conversion.factor = spot.size / scalefactors$spot_diameter_fullres
  spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size /
                                 2)
  
  # Create cell chat object
  message(paste('Creating CellChat object for sample:', samples[x]))
  cellchat <-
    createCellChat(
      object = data.input,
      meta = meta,
      group.by = "labels",
      datatype = "spatial",
      coordinates = data.matrix(spatial.locs),
      spatial.factors = spatial.factors
    )
  
  # Set database
  load('cellphonedb/CellChatDB.human_user.rda')
  cellchat@DB <- db.new
  
  # Subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat)
  future::plan("multisession", workers = 4)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
  
  # Compute communication probabilities
  message(paste('Computing communication probabilities for sample:', samples[x]))
  cellchat <-
    computeCommunProb(
      cellchat,
      type = "truncatedMean",
      trim = 0.1,
      distance.use = TRUE,
      interaction.range = 250,
      scale.distance = 0.01,
      contact.dependent = TRUE,
      contact.range = 100
    )
  
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # Infer communcation at pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Calculate aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  
  # Save cellchat object
  message(paste('Saving cellchat object for sample:', samples[x]))
  saveRDS(cellchat, file = paste0(sampledir, '/', samples[x], '_cellchat.rds'))
  
}

# Load the above generated cellchat objects 
rdsobjs <- as.list(list.files(path = 'cellChat_outs', 
                              pattern = '*.rds',
                              recursive = T,
                              full.names = T))
rdslist <- lapply(rdsobjs, function(x){
  readRDS(x)
})
names(rdslist) <- gsub(
  '_cellchat.rds',
  '',
  basename(list.files(
    path = 'cellChat_outs',
    pattern = '*.rds',
    recursive = T,
    full.names = T
  ))
)

# Set outdir
outdir <- 'cellChat_outs/'

# Visualise interactions for all samples in loop ----
for (ITER in names(rdslist)) {
  
  # Set sample outdir
  sampledir <- file.path(outdir, ITER)
  if (!dir.exists(sampledir)) {
    dir.create(sampledir)
  }
  
  if (sum(rowSums(rdslist[[ITER]]@net$count)) == 0) {
    message(paste0('No interactions for ', ITER, ': SKIPPING...'))
    next
  } 
  
  # Plot number of interactions as circle plot
  groupSize <- as.numeric(table(rdslist[[ITER]]@idents))
  pdf(
    file = paste0(sampledir, '/', ITER, '_num_intrns.pdf'),
    height = 4,
    width = 4
  )
  netVisual_circle(
    rdslist[[ITER]]@net$count,
    vertex.weight = rowSums(rdslist[[ITER]]@net$count),
    weight.scale = T,
    label.edge = F,
    title.name = NULL
  )
  dev.off()
  
  # Plot strength of interactions as circle plot
  pdf(
    file = paste0(sampledir, '/', ITER, '_strgth_intrns.pdf'),
    height = 4,
    width = 4
  )
  netVisual_circle(
    rdslist[[ITER]]@net$weight,
    vertex.weight = rowSums(rdslist[[ITER]]@net$weight),
    weight.scale = T,
    label.edge = F,
    title.name = NULL
  )
  dev.off()
  
  # Plot number and strength of interactions as heatmaps
  pdf(
    file = paste0(sampledir, '/', ITER, '_num_intrns_heat.pdf'),
    height = 4,
    width = 4
  )
  print(netVisual_heatmap(rdslist[[ITER]],
                    measure = "count",
                    color.heatmap = "Blues"))
  dev.off()
  
  pdf(
    file = paste0(sampledir, '/', ITER, '_strgth_intrns_heat.pdf'),
    height = 4,
    width = 4
  )
  print(netVisual_heatmap(rdslist[[ITER]],
                    measure = "weight",
                    color.heatmap = "Blues"))
  dev.off()
  
  
  # All pathways showing significant interaction:
  pathways.show.all <- rdslist[[ITER]]@netP$pathways
  write.csv(pathways.show.all, file = paste0(sampledir, '/', ITER, '_pathways.csv'))
  
  # Export the database for reference
  #write.csv(rdslist[[1]]@DB[["interaction"]], 'cellChatDB.csv')
  
  # It is possible to loop through the above and generate plots for all
  # Define pathways of interest
  # poi <- c('TAC',
  #          'OPIOID',
  #          'GALANIN',
  #          'APELIN',
  #          'IL6',
  #          'CXCL',
  #          'CCL',
  #          'IL1',
  #          'NPR1',
  #          'CALCR')
  # 
  # # Plot the interactions between clusters as circle plot for specific pathways
  # for (x in poi) {
  #   if (x %in% pathways.show.all) {
  #     pdf(
  #       file = paste0(ITER, '_', x, '_circle.pdf'),
  #       height = 4,
  #       width = 4
  #     )
  #     netVisual_aggregate(
  #       rdslist[[ITER]],
  #       signaling = x,
  #       layout = "circle",
  #       color.use = NULL,
  #       sources.use = NULL,
  #       targets.use = NULL,
  #       idents.use = NULL
  #     )
  #     dev.off()
  #     
  #     # Plot the interactions between clusters as heatmap for specific pathways
  #     pdf(
  #       file = paste0(ITER, '_', x, '_heat.pdf'),
  #       height = 4,
  #       width = 4
  #     )
  #     print(netVisual_heatmap(rdslist[[ITER]], signaling = x, color.heatmap = "Reds"))
  #     dev.off()
  #   } else {
  #     next
  #   }
  # }
  # 
  # # Plot all the significant interactions (L-R pairs) from some cell groups
  # # Define cell groups by 'sources.use' to other cell groups (defined by 'targets.use')
  # # Put cluster number in quotes to get exact match
  # 
  # # Check idents
  # sampleIdents <- unique(rdslist[[ITER]]@idents)
  # 
  # # Filter poi to only valid pathways
  # poi <- poi[poi %in% pathways.show.all]
  # 
  # # Define sources and targets - source = 10 - targets all others
  # sources <- '10'
  # targets <-
  #   c(
  #     '1',
  #     '2',
  #     '3',
  #     '4',
  #     '5',
  #     '6',
  #     '7',
  #     '8',
  #     '9',
  #     '10',
  #     '11',
  #     '12',
  #     '13',
  #     '14',
  #     '15',
  #     '16',
  #     '17',
  #     '18',
  #     '19'
  #   )
  # # Bubble plot
  # pdf(
  #   file = paste0(ITER, '_10source_bubble.pdf'),
  #   height = 4,
  #   width = 4
  # )
  # print(netVisual_bubble(
  #   rdslist[[ITER]],
  #   signaling = poi,
  #   sources.use = sources,
  #   targets.use = targets,
  #   remove.isolate = FALSE,
  #   font.size = 4
  # ))
  # dev.off()
  # 
  # 
  # # Define sources and targets - target = 10 - sources all others
  # targets <- '10'
  # sources <-
  #   c(
  #     '1',
  #     '2',
  #     '3',
  #     '4',
  #     '5',
  #     '6',
  #     '7',
  #     '8',
  #     '9',
  #     '10',
  #     '11',
  #     '12',
  #     '13',
  #     '14',
  #     '15',
  #     '16',
  #     '17',
  #     '18',
  #     '19'
  #   )
  # # Bubble plot
  # pdf(
  #   file = paste0(ITER, '_10target_bubble.pdf'),
  #   height = 4,
  #   width = 4
  # )
  # print(netVisual_bubble(
  #   rdslist[[ITER]],
  #   signaling = poi,
  #   sources.use = sources,
  #   targets.use = targets,
  #   remove.isolate = FALSE,
  #   font.size = 4
  # ))
  # dev.off()
  # 
  # # Define sources and targets - source = 2 - targets all others
  # sources <- '2'
  # targets <-
  #   c(
  #     '1',
  #     '2',
  #     '3',
  #     '4',
  #     '5',
  #     '6',
  #     '7',
  #     '8',
  #     '9',
  #     '10',
  #     '11',
  #     '12',
  #     '13',
  #     '14',
  #     '15',
  #     '16',
  #     '17',
  #     '18',
  #     '19'
  #   )
  # # Bubble plot
  # pdf(
  #   file = paste0(ITER, '_2source_bubble.pdf'),
  #   height = 4,
  #   width = 4
  # )
  # print(netVisual_bubble(
  #   rdslist[[ITER]],
  #   signaling = poi,
  #   sources.use = sources,
  #   targets.use = targets,
  #   remove.isolate = FALSE,
  #   font.size = 4
  # ))
  # dev.off()
  # 
  # # Define sources and targets - target = 2 - sources all others
  # targets <- '2'
  # sources <-
  #   c(
  #     '1',
  #     '2',
  #     '3',
  #     '4',
  #     '5',
  #     '6',
  #     '7',
  #     '8',
  #     '9',
  #     '10',
  #     '11',
  #     '12',
  #     '13',
  #     '14',
  #     '15',
  #     '16',
  #     '17',
  #     '18',
  #     '19'
  #   )
  # # Bubble plot
  # pdf(
  #   file = paste0(ITER, '_2target_bubble.pdf'),
  #   height = 4,
  #   width = 4
  # )
  # print(netVisual_bubble(
  #   rdslist[[ITER]],
  #   signaling = poi,
  #   sources.use = sources,
  #   targets.use = targets,
  #   remove.isolate = FALSE,
  #   font.size = 4
  # ))
  # dev.off()
  
}










