# User edit here:

# Where you want outputs to be saved 
outdir <- '/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/NanostringGeoMX/Projects/Fanti_Silvia/GC-SF-10611/Analysis/Run_R_Script_Remote_Test'

# Folder containing DCC files output from NGS pipeline
dccpath <- '/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/NanostringGeoMX/Projects/Fanti_Silvia/GC-SF-10611/Analysis/NGS_Pipeline_Outs_v2'

# Folder containing pkc file
pkcpath <- '/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/NanostringGeoMX/Projects/Fanti_Silvia/GC-SF-10611/Analysis/R_Analysis_v2'

# SampleAnnotationFile: This file needs to be an xlsx worksheet
# Needs to be made from the LabWorksheet.txt file downloaded with the config.ini file from machine
# Delete the top rows of the file so only phenotypic data is left with colnames
# Call the sheet Template
SampleAnnotationFile <- '/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/NanostringGeoMX/Projects/Fanti_Silvia/GC-SF-10611/Analysis/R_Analysis_v2/GC-SF-10611_20230901T1526_LabWorksheetMOD.xlsx'

# Names of the aoi and roi columns in the above spreadsheet
protColNames <- c("aoi", "roi")

# Name of the column with the panel information
experimentColNames <- 'panel'

# Do not edit beyond here

# Set R library path
.libPaths('/data/WHRI-GenomeCentre/software/R-library-GeoMx')

# Load libraries
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
  path = dccpath,
  pattern = "*.dcc",
  full.names = TRUE,
  recursive = TRUE
)

# PKC file downloaded from nanostring website
# These contain probe assay meta data - download from web
PKCFiles <- dir(
  path = pkcpath,
  pattern = ".pkc",
  full.names = TRUE,
  recursive = TRUE
)

# Load data
fullData <- readNanoStringGeoMxSet(
  dccFiles = DCCFiles,
  pkcFiles = PKCFiles,
  phenoDataFile = SampleAnnotationFile,
  phenoDataSheet = "Template",
  phenoDataDccColName = "Sample_ID",
  protocolDataColNames = protColNames,
  experimentDataColNames = experimentColNames
)

# STUDY DESIGN ----

# Modules used 
pkcs <- annotation(fullData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

# QC & PRE-PROCESSING ----

# Shift counts to one - to enable downstream transformations
fullData <- shiftCountsOne(fullData, useDALogic = TRUE)

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

# Calculate the negative probe geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(fullData), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 

colnames(negativeGeoMeans) <- 'NTC_Neg_Geo_Means'

protocolData(fullData)[["NegGeoMean"]] <- negativeGeoMeans

# Explicitly copy the negative probe geometric means from sData to pData
pData(fullData)[, "NegGeoMean"] <- sData(fullData)[["NegGeoMean"]]

# Detach negative probe geometric mean columns ahead of aggregateCounts call
pData(fullData) <- pData(fullData)[, !colnames(pData(fullData)) %in% "NegGeoMean"]

# Show all NTC values, Freq = # of Segments with a given NTC count:
# NTC = NO TEMPLATE CONTROL - establishes the level at which counts in the NTC will be flagged
# NTC is used to detect contamination in the library prep
# Save with other metrics read_summary.csv
summaryTable <- data.frame(
  sample_ID = fullData@protocolData@data$SampleID,
  raw_Reads = fullData@protocolData@data$Raw,
  trimmed_Reads = fullData@protocolData@data$Trimmed,
  stiched_Reads = fullData@protocolData@data$Stitched,
  aligned_Reads = fullData@protocolData@data$Aligned,
  deduplicated_Reads = fullData@protocolData@data$DeduplicatedReads,
  NTC_GeoMean = fullData@protocolData@data$NegGeoMean
)
write.csv(summaryTable, 'read_summary.csv')

# Save segment QC summary
write.csv(QC_Summary, "segment_qc_summary.csv")

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

# Save raw probe count matrix
write.csv(exprs(fullData), 'Raw_Probe_Counts.csv')

# CREATE GENE LEVEL COUNTS ----

# Collapse to targets
target_fullData <- aggregateCounts(fullData)

# Export raw counts (shifted to 1)
raw_counts <- exprs(target_fullData)
write.csv(raw_counts, 'Raw_Gene_Counts.csv')