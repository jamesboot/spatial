#!/bin/bash 

## Running Spaceranger mkfastq to demultiplex raw sequence data (.bcl) into fastq files for downstream analysis. 

##### Set queue options
#$ -M j.boot@qmul.ac.uk
#$ -m bes
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=40:0:0
#$ -l h_vmem=10G
#$ -N spaceranger_mkfastq
#$ -j y
#####

### load the module ## 
module load spaceranger

### load bcl2fastq2 ##
module load bcl2fastq2

## run Spaceranger mkfastq ##
## --id is the name of the folder that will be created by Spaceranger mkfastq ##
## --run is the full path to the sequencing run folder ##
## --csv is the full path to the .csv file containing sample and index information
spaceranger mkfastq --id=GC-DS-9511_fastq_full \
                   --run=/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/10x_Chromium/Schalke_Diana/GC-DS-9511/Data/211102_NS500784_0827_AHHVM5BGXK \
                   --csv=/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/10x_Chromium/Schalke_Diana/GC-DS-9511/Analysis/GC-DS-9511-mkfastq-simple.csv
