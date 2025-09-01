# GeoMx_Analysis
- All scripts associated with GeoMX DSP data analysis.
- See OneNote SOP for step by step details.

## 1. rename-fastqs.sh
- FASTQ files must be in the correct name format, or the NGS pipeline script will fail.
- Required format: <SAMPLE_ID>_Sx_L001_R[1,2]_001.fastq.gz
- Depending on the sequencer, L001 is sometimes missing - rename-fastqs.sh is designed to fix this 

## 2. GeoMx-NGS-Pipeline.sh
- Script for pre-processing of raw data from FASTQ files to DCC files.
- DCC files can be imported onto the DSP or into R for further analysis.

## 3. GeoMx_countMat.R & GeoMx_CountMat_sub.sh
- GeoMx_countMat.R is an R script which will generate:
    - read_summary.csv - CSV file containing QC metrics per ROI
    - segment_qc_summary.csv - CSV file containing summary of how many segments (ROIs) passed or failed different QC thresholds
        - QC thresholds are in the script - defaults used
    - Raw_Probe_Counts.csv - count matrix at probe level 
    - Raw_Gene_Counts.csv - count matrix at gene level
- GeoMx_countMat.R uses an R environment set up in the WHRI-GenomeCentre/softwares directory - it is specified in the script.
- GeoMx_countMat_sub.sh is a submission script to run the R script on Apocrita 

## Other:

### GeoMx_Analysis.R
- R script written by James Boot for analysis of GeoMx DCC files. 
- This script will generate or perform, QA/QC, filtering, raw count matrices, normalised count matrices. 
- PKC file needed for this script. See PKC files in this repository.
- This script is designed to run in an interative R session in R studio. Do not run on Apocrita.

### GeoMx_Manuals
- Folder contains various manuals from Nanostring for DSP operation and data analysis.

### Hsa_WTA_v1.0.pkc
- Human whole transcriptome PKC file needed for analysis with GeoMx_Analysis.R script.
