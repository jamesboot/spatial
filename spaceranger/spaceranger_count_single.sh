#!/bin/bash 

## Running spaceranger for count analysis

#### Set queue options
#$ -M e.wozniak@qmul.ac.uk			# Add your email address
#$ -m bes					# Request email notifications
#$ -cwd						# Set to current working directory
#$ -pe smp 4					# Request xx cores
#$ -l h_rt=24:0:0				# Request xx running time
#$ -l h_vmem=8G					# Request xx memory
#$ -N count-AB-9124				# Name the job
#$ -j y						# Join the standard error/output stream
####

## load the module
module load spaceranger
    
## Commands to run count 
# --id is the name of the output folder created by spaceranger count
# --fastqs is the directory containing the fastq files produced by spaceranger mkfastq
# --transcriptome is the reference to align the data to
# --sample is the sample ID of the fastq files
# --image is the directory of the image file (TIFF)

spaceranger count --id=GC-AB-9124-GEX_full \
                 --fastqs=/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/10x_Chromium/Biddle_Adrian/GC-AB-9124/Analysis/GC-AB-9124_fastq_full/outs/fastq_path/HGJY5BGXH \
                 --transcriptome=/data/WHRI-GenomeCentre/Genome/cellranger_5.0_references/human/refdata-gex-GRCh38-2020-A/ \
                 --sample=GC-AB-9124-GEX
                 --image=/path/to/image/file.tif
                 --slide=<<SLIDEREF>>
                 --area=<<SLIDE AREA>>
