setwd('/Users/jamesboot/Documents/9.Genome Centre Files/GC-WG-10782/')

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(ggplot2)
library(knitr)
library(scales)
library(reshape2) 
library(cowplot)
library(umap)
library(Rtsne)
library(pheatmap)
library(ggbreak)
library(dplyr)
library(magrittr)
library(kableExtra)

# DATA LOADING ----

# Automatically list files in each directory for use
# DCC file made from the NGS pipeline
# These contain count data and sequencing quality meta data
DCCFiles <- dir(
  path = '/Users/jamesboot/Documents/9.Genome Centre Files/GC-WG-10782/NTC_data/GeoMx_NGS_Pipeline_Outs/',
  pattern = "*.dcc",
  full.names = TRUE,
  recursive = TRUE
)

# PKC file downloaded from nanostring website
# These contain probe assay meta data - download from web
PKCFiles <- dir(
  path = '/Users/jamesboot/Documents/9.Genome Centre Files/GC-WG-10782',
  pattern = ".pkc",
  full.names = TRUE,
  recursive = TRUE
)

# This file needs to be an xlsx worksheet
# Needs to be made from the LabWorksheet.txt file downloaded with the config.ini file from machine
# Delete the top rows of the file so only phenotypic data is left with colnames
# Call the sheet Template
SampleAnnotationFile <- 'NTC_data/GC-WG-10782_Negative_Collection_20240409T1335_LabWorksheet.xlsx'

# Would you like to investigate the no template control counts?
NTC <- FALSE

if (NTC == TRUE) {
  # Load data
  fullData <- readNanoStringGeoMxSet(
    dccFiles = DCCFiles,
    pkcFiles = PKCFiles,
    phenoDataFile = SampleAnnotationFileMOD,
    phenoDataSheet = "Template",
    phenoDataDccColName = "Sample_ID",
    protocolDataColNames = c("aoi", "roi"),
    experimentDataColNames = c("panel")
  )
  
  # NO TEMPLATE CONTROL INVESTIGATION ----
  
  # Extract counts
  raw_counts <- fullData@assayData[["exprs"]]
  
  sum(raw_counts[, c('DSP-1001660010992-B-A01.dcc')])
  sum(raw_counts[, c('DSP-1001660010993-C-A01.dcc')])
  
  write.csv(raw_counts, file = 'raw_probe_counts_inclNTC.csv')
  
  summaryTable <- data.frame(
    sample_ID = fullData@protocolData@data$SampleID,
    raw_Reads = fullData@protocolData@data$Raw,
    trimmed_Reads = fullData@protocolData@data$Trimmed,
    stiched_Reads = fullData@protocolData@data$Stitched,
    aligned_Reads = fullData@protocolData@data$Aligned,
    deduplicated_Reads = fullData@protocolData@data$DeduplicatedReads
  )
  write.csv(summaryTable, 'read_summary.csv')
  
  # Heatmap of counts
  library(gplots)
  tiff(
    filename = 'NTC_Heatmap.tiff',
    width = 9,
    height = 9,
    res = 300,
    units = 'in'
  )
  heatmap.2(
    raw_counts,
    dendrogram = 'col',
    col = 'bluered',
    scale = 'row',
    trace = 'none',
    labRow = F,
    cexCol = 0.25
  )
  dev.off()
  
} else {
  
  # Load data
  fullData <- readNanoStringGeoMxSet(
    dccFiles = DCCFiles,
    pkcFiles = PKCFiles,
    phenoDataFile = SampleAnnotationFile,
    phenoDataSheet = "Template",
    phenoDataDccColName = "Sample_ID",
    protocolDataColNames = c("Aoi", "Roi"),
    experimentDataColNames = c("Panel")
  )
  
}

# STUDY DESIGN ----

# Modules used 
pkcs <- annotation(fullData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

# QC & PRE-PROCESSING ----

# Shift counts to one - to enable downstream transformations
fullData <- shiftCountsOne(fullData, useDALogic = TRUE)

# Extract the de-duplicated read counts and other metrics
summaryTable <- data.frame(
  sample_ID = fullData@protocolData@data$SampleID,
  raw_Reads = fullData@protocolData@data$Raw,
  trimmed_Reads = fullData@protocolData@data$Trimmed,
  stiched_Reads = fullData@protocolData@data$Stitched,
  aligned_Reads = fullData@protocolData@data$Aligned,
  deduplicated_Reads = fullData@protocolData@data$DeduplicatedReads
)
write.csv(summaryTable, 'read_summary.csv')

# SEGMENT QC ----

# Default QC cutoffs are commented in () 
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 80,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 10,   # Minimum negative control counts (10)
       maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 100,         # Minimum # of nuclei estimated (100)
       minArea = 5000)         # Minimum segment area (5000)

fullData <-
  setSegmentQCFlags(fullData, 
                    qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(fullData)[["QCFlags"]]

flag_columns <- colnames(QCResults)

QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))

QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})

QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

# Change here what we want to separate the data by - use the colnames of the LabWorksheet!
col_by <- "Segment"

# Visualise QC metrics 
# Graphical summaries of segment QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

QC_histogram(sData(fullData), "Trimmed (%)", col_by, 80)
QC_histogram(sData(fullData), "Stitched (%)", col_by, 80)
QC_histogram(sData(fullData), "Aligned (%)", col_by, 75)
QC_histogram(sData(fullData), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
QC_histogram(sData(fullData), "Area", col_by, 1000, scale_trans = "log10")
# Need to add nuclei else where - needs to be exported from DSP - future work to add this to the script
#QC_histogram(sData(fullData), "nuclei", col_by, 20)

# Calculate the negative probe geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(fullData), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 

protocolData(fullData)[["NegGeoMean"]] <- negativeGeoMeans

# Explicitly copy the negative probe geometric means from sData to pData
negCols <- paste0("NegGeoMean_", modules)

pData(fullData)[, negCols] <- sData(fullData)[["NegGeoMean"]]

for(ann in negCols) {
  plt <- QC_histogram(pData(fullData), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}

# Detach negative probe geometric mean columns ahead of aggregateCounts call
pData(fullData) <- pData(fullData)[, !colnames(pData(fullData)) %in% negCols]

# Show all NTC values, Freq = # of Segments with a given NTC count:
# NTC = NO TEMPLATE CONTROL - establishes the level at which counts in the NTC will be flagged
# NTC is used to detect contamination in the library prep
kable(table(NTC_Count = sData(fullData)$NTC),
      col.names = c("NTC Count", "# of Segments"))

kable(QC_Summary, caption = "QC Summary Table for each Segment") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  cat(., file = "segment_qc_summary.html")

# Removed flagged data - we do not do this!
#fullData <- fullData[, QCResults$QCStatus == "PASS"]

# Subsetting our dataset has removed samples which did not pass QC
dim(fullData)

# PROBE QC ----

# Generally keep the qcCutoffs parameters unchanged. 
# Set removeLocalOutliers to FALSE if you do not want to remove local outliers
fullData <- setBioProbeQCFlags(fullData, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = FALSE)

ProbeQCResults <- fData(fullData)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))

# Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(fullData, 
         fData(fullData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(fullData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)

fullData <- ProbeQCPassed 

# Check how many unique targets the object has
length(unique(featureData(fullData)[["TargetName"]]))

write.csv(exprs(fullData), 'raw_probe_count_mat.csv')

# CREATE GENE LEVEL COUNTS ----

# Collapse to targets
target_fullData <- aggregateCounts(fullData)
dim(target_fullData)

# Visualise first 5 genes in first 10 samples as table
exprs(target_fullData)[1:5, 1:10]

# Export raw counts (shifted to 1)
raw_counts <- exprs(target_fullData)
write.csv(raw_counts, 'GC-WG-10782_Raw_Counts.csv')

# Normalise here prior to further filtering
norm_fullData <- normalize(target_fullData ,
                           norm_method = "quant", 
                           desiredQuantile = .75,
                           toElt = "q_norm")

norm_counts1 <- norm_fullData@assayData[["q_norm"]]
write.csv(norm_counts1, 'GC-WG-10782_Q3norm_Counts.csv')

# LIMIT OF QUANTIFICATION ----

# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_fullData))

for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_fullData)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_fullData)[, vars[1]] * 
             pData(target_fullData)[, vars[2]] ^ cutoff)
  }
}

pData(target_fullData)$LOQ <- LOQ

# FILTERING ---- 

# Filter out segments and genes with low signal
LOQ_Mat <- c()

for(module in modules) {
  ind <- fData(target_fullData)$Module == module
  Mat_i <- t(esApply(target_fullData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}

# Ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_fullData)$TargetName, ]

# Low signal segments will have a small fraction of panel genes detected
# Above the LOQ - visualise distribution of segments relative to detected genes 
# Save detection rate information to pheno data
pData(target_fullData)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)

pData(target_fullData)$GeneDetectionRate <-
  pData(target_fullData)$GenesDetected / nrow(target_fullData)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_fullData)$DetectionThreshold <- 
  cut(pData(target_fullData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# Stacked bar plot of different cut off points (1%, 5%, 10%, 15%)
# This is the percentage of genes detected and how many segments meet that threshold
ggplot(pData(target_fullData),
       aes(x = DetectionThreshold)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")

# Show above as a table as well
kable(table(pData(target_fullData)$DetectionThreshold))

# Generally they cut out anything less than 10% - here we use 5% to be lenient
target_fullData <-
  target_fullData[, pData(target_fullData)$GeneDetectionRate >= .05]

dim(target_fullData)

# Calculate gene detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_fullData)]

fData(target_fullData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)

fData(target_fullData)$DetectionRate <-
  fData(target_fullData)$DetectedSegments / nrow(pData(target_fullData))

# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))

plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_fullData)$DetectionRate >= x)}))

plot_detect$Rate <- plot_detect$Number / nrow(fData(target_fullData))

rownames(plot_detect) <- plot_detect$Freq

# How many genes are detect in certain percentages of segments?
ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

# Subset to target genes detected in at least 5% of the samples.
# Normal cut off is anywhere from 5-20%
# Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_fullData), CodeClass == "Negative")

neg_probes <- unique(negativeProbefData$TargetName)

target_fullData <- 
  target_fullData[fData(target_fullData)$DetectionRate >= 0.05 |
                    fData(target_fullData)$TargetName %in% neg_probes, ]

dim(target_fullData)

# Export the filtered raw counts
write.csv(exprs(target_fullData),
          'GC-MS-10537_Filtered_Raw_Counts.csv')

# DATA NORMALISATION ----

# Q3 norm is recommended for NGS readouts WTA/CTA with or without custom spike-ins
# Desired quartile is 75th
# Q3 divides the counts in one segment by the 3rd quartile value for that segment...
# Then multiplies that value by the geometric mean of the 3rd quartile value of all segments
# Perform after filtering
norm_filtData <- normalize(target_fullData ,
                           norm_method = "quant", 
                           desiredQuantile = .75,
                           toElt = "q_norm")

# Export the normalised counts
norm_counts2 <- norm_filtData@assayData[["q_norm"]]
write.csv(norm_counts2, 'GC-MS-10537_Filtered_Q3Norm_Counts.csv')
