# Script for running NoButter package on CosMx data
# Analyses read distribution over z-stacks
# Removes noisey z-stacks

# Snapshot env
renv::snapshot()

# Load packages
library(NoButter)
library(ggplot2)
library(dplyr)
library(grid)
library(patchwork)
library(reshape2)
library(pheatmap)
library(data.table)

# Do we neeed to combined raw transcripts from each FOV?
# Only needs ot be done once, takes long time
combineTx <- F
if (combineTx == T) {
  # Make a list of the 2 slide locations to iterate through
  locs <- c(
    "data/AtoMx_Oct24/RawFiles/GCIN10692_TMA2/AnalysisResults/t8edqe0wnx/",
    "data/AtoMx_Oct24/RawFiles/GCIN10692_TMA3/AnalysisResults/di3i9kvpg6/"
  )
  sample <- c("GCIN10692_TMA2", "GCIN10692_TMA3")
  
  # Applying function to combine raw transcripts data from each FOV into one .csv file
  for (x in 1:length(locs)) {
    CombineRawTranscripts(
      in_dir = locs[x],
      out_dir = locs[x],
      file_name = paste0(sample[x], "_combined_transcripts.csv")
    )
  }
}

# Define parameters & inputs depending which slide is being analysed
slideID <- 'GCIN10692_TMA2'
slideNo <- 1
txCsv <- "data/AtoMx_Oct24/RawFiles/GCIN10692_TMA2/AnalysisResults/t8edqe0wnx/GCIN10692_TMA2_combined_transcripts.csv"
metaCsv <- "data/AtoMx_Oct24/flatFiles/GCIN10692_TMA2/GCIN10692_TMA2_metadata_file.csv.gz"
polyCsv <- "data/AtoMx_Oct24/flatFiles/GCIN10692_TMA2/GCIN10692_TMA2-polygons.csv.gz"
flatFileOuts <- "data/AtoMx_Oct24/sample_dir_formatted/GCIN10692_TMA2/"
outDir <- paste0('NoButter_Outs/', slideID)

# Define output directory for plots
if (!dir.exists(outDir)){
  dir.create(outDir)
}

# view transcripts file
transcripts <- read.csv(txCsv) 
head(transcripts)

# file dimensions
dim(transcripts)

# view what information is contained within the columns of the transcripts file
colnames(transcripts)

# removing redundant column 
transcripts$X <- NULL

# add slide id column (breaks later without it)
transcripts$Slide = slideNo

# check how many z stacks present in the dataset
unique(transcripts$z)

# Plot Transcript Counts

# plot number of transcripts per FOV - line at 30,000
pdf(file = paste0(outDir, '/', slideID, '_Tx_perFOV.pdf'))
p <- PlotTranscriptCounts(transcripts)
plot(p)
dev.off()

# code class distribution
pdf(file = paste0(outDir, '/', slideID, '_Tx_class.pdf'))
p <- PlotTranscriptByCodeClass(transcripts)
plot(p)
dev.off()

# Z-Stack investigation 
# Number of transcripts per z stack

# investigating number of transcripts per z stack
zstack_count_tbl <- table(transcripts$z)
zstack_count_tbl

# visualizing transcripts per z stack
counts <- setNames(as.data.frame(zstack_count_tbl), c("z", "Count"))
pdf(file = paste0(outDir, '/', slideID, '_Tx_perZ.pdf'))
p <- ggplot(data = counts, aes(x = z, y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Transcripts per Z Stack", x = "Z Stack", y = "Count") +
  scale_y_continuous(labels = scales::comma) +
  theme_custom()
plot(p)
dev.off()

# proportion of transcripts per z stack
zstack_prop_tbl <- round(table(transcripts$z)*100/nrow(transcripts), 3)
zstack_prop_tbl

# proportion of transcripts per z stack
prop <- setNames(
  as.data.frame(round(table(transcripts$z) * 100 / nrow(transcripts), 3)),
  c("z_stack", "Proportion")
)

pdf(file = paste0(outDir, '/', slideID, '_Tx_perZprop.pdf'))
p <- ggplot(prop, aes(x = 2, y = Proportion, fill = z_stack)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Proportion of transcripts per z stack") +
  theme_void() +
  theme(legend.position = "right")
plot(p)
dev.off()

# Plot distribution as heatmap
VisualizeZStackDistribution(transcripts, fig_save = T, out_dir = outDir)
dev.off()

# Plot distribution as heatmap in binary
VisualizeZStackDistribution(
  transcripts,
  binarised = TRUE,
  fig_save = T,
  out_dir = outDir
)
dev.off()

# Specify the fovs to keep
keep_zstack <- 0:7

# proportion of endogenous transcripts kept
PlotTranscriptLossProportion(transcripts, 
                             keep_z_stacks = keep_zstack,
                             fig_save = T,
                             out_dir = outDir)
dev.off()

# filter and summarise data
summary_data <- transcripts %>%
  mutate(Category = ifelse(z %in% keep_zstack, "Kept", "Lost")) %>%
  group_by(Category) %>%
  summarise(
    Count = n(),
    Percentage = (n() / nrow(transcripts)) * 100
  ) %>%
  ungroup()
print(summary_data)

# Here we look at both the impact of transcript removal 
# alongside transcripts that were not assigned to a cell (cellID=0)
# In an example FOV
exampleFOV <- 57

p <- transcripts %>% 
  filter(fov == exampleFOV) %>%
  ggplot(aes(x=x, y=y, colour=CellId>0)) + 
  geom_point(size=0.1)  + 
  facet_wrap(~z %in% keep_zstack)
ggsave(plot = p, filename = paste0(outDir, '/', slideID, '_TxRemoval1.pdf'))

# View distribution transcripts that will be kept in specific
p <- transcripts %>% 
  filter(fov == exampleFOV) %>%
  ggplot(aes(x=x, y=y, z %in% keep_zstack)) + 
  geom_point(size=0.1)  + 
  facet_wrap(~z %in% keep_zstack)
ggsave(plot = p, filename = paste0(outDir, '/', slideID, '_TxRemoval2.pdf'))

# Another view - in cells vs not in cells
p <- transcripts %>% filter(fov == exampleFOV) %>%
  ggplot(aes(x = x, y = y, colour = CellId > 0)) + 
  geom_point(size = 0.1)  + 
  facet_wrap( ~ CellId > 0)
ggsave(plot = p, filename = paste0(outDir, '/', slideID, '_TxRemoval3.pdf'))

# Updating flat files
# Clean transcripts file
# Check for presence of transcripts that fall outside of cell boundaries (cellID==0)
table(transcripts$CellId == 0)

# Calculate the count of transcripts within each specified z stack range
transcript_counts <- transcripts %>%
  mutate(z_stacks = case_when(
    z %in% keep_zstack ~ "To keep",
    (!z %in% keep_zstack ~ "To remove")
  )) %>%
  group_by(z_stacks) %>%
  summarise(count = n())
transcript_counts

# Remove transcripts from z stacks 8 and 9 and remove transcripts not contained within cell boundaries
transcripts_zclean <- transcripts %>%
  filter(z %in% keep_zstack & CellId != 0)
dim(transcripts_zclean)
head(transcripts_zclean)

# saving cleaned Tx data
write.csv(transcripts_zclean,
          paste0(flatFileOuts, '/', slideID, "_tx_zclean.csv"))

# Generate clean gene expression matrix
# PrepareGEXMatrix() function key steps: 
# - remove cells with id == 0 (not within segmentation mask) 
# - remove SystemControl probes 
# - keep only z stack between x:y 
# - create unique cell id "c_slide_fov_cellID" - crosstab between cell id and target

# Using PrepareGEXMatrix function to generate a clean gene expression matrix
GEX_clean_z_stack <- PrepareGEXMatrix(transcripts, 
                                      z_stacks = keep_zstack, 
                                      save_matrix = F)

# Clean metadata and polygons 

# The *UpdateMetadataAndPolygons()* function is designed to update and synchronize metadata and polygon files with...
# cleaned gene expression (GEX) data, primarily by filtering and reordering them to match the cleaned dataset.
# Key steps in function: 
# - Aligns the metadata and polygons with the cleaned GEX data, including only those entries corresponding to cells present in the GEX data and arranging them in the same order
# - Recalculates important metrics for the metadata based on the cleaned GEX data, such as total RNA counts, number of features (genes) detected, and counts for both regular and negative control probes.
metadata <- fread(metaCsv,)
polygons <- fread(polyCsv,)

UpdateMetadataAndPolygons(clean_GEX = GEX_clean_z_stack,
                          metadata = metadata,
                          polygons = polygons,
                          save_flat_files = TRUE,
                          output_dir = flatFileOuts)

# Add column called cell_ID for importing into squidpy - should just be numeric
strs <- strsplit(rownames(GEX_clean_z_stack), "_")
cells <- unlist(lapply(strs, function(x){
  x[4]
}))
GEX_clean_z_stack$cell_ID <- as.numeric(cells)

# Add column for fov fr importing into squidpy - should be just numeric
strs <- strsplit(rownames(GEX_clean_z_stack), "_")
fovs <- unlist(lapply(strs, function(x){
  x[3]
}))
GEX_clean_z_stack$fov <- as.numeric(fovs)

# Remove rownames
rownames(GEX_clean_z_stack) <- NULL

# Save matrix
write.csv(GEX_clean_z_stack,
          file = paste0(flatFileOuts, 
                        '/',
                        slideID,
                        '_exprMat_zclean.csv'),
          quote = F,
          row.names = F)

test1 <- read.csv('../../data/AtoMx_Oct24/sample_dir_formatted/GCIN10692_TMA2/GCIN10692_TMA2_exprMat_file.csv')
test2 <- read.csv('../../data/AtoMx_Oct24/sample_dir_formatted/GCIN10692_TMA2/GCIN10692_TMA2_exprMat_zclean.csv')
