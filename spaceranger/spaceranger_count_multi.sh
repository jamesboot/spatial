#!/bin/bash 

#### SECTION 1 - Set queue options

#$ -N spaceranger-count           # Name the job
#$ -M j.boot@qmul.ac.uk           # Add your email
#$ -m bes                         # Request email notifications
#$ -cwd                           # Set to current working directory
#$ -j y                           # Join standard error / log stream
#$ -pe smp 16                     # Request 16 cores
#$ -l h_vmem=8G                   # Request 8G RAM
#$ -l h_rt=48:00:00               # Request 48 hours run time - each sample takes anywhere from 4-11 hours to run depending on size (6-7 seems average)
#$ -t 1-4                         # Set -t to number of array tasks to run - change this based on number of samples
#$ -l highmem					  # Join high memory queue



#### SECTION 2 - Load the module

module load spaceranger

#### SECTION 3 - Specify the parameters file and arguements to be used later
# parameters.txt is a text file which has mulitple lines - each sample is on a different line - needs to be prepared by user
# On each iteration of the array the script will use a different line i.e. iteration 1 = line 1, iteration 2 = line 2 etc.
# Each line contains multiple pieces of information for each sample separated by a colon
# First part of the line is the sample name (SAMPLE), then a colon, then the list of fastq directories (FASTQ_DIR) seperated by a comma

INPUT_ARGS=$(sed -n "${SGE_TASK_ID}p" GC-DS-9511-parameters.txt)
SAMPLE=$(echo $INPUT_ARGS | cut -d : -f 1)
FASTQ_DIR=$(echo $INPUT_ARGS | cut -d : -f 2)
IMAGE_DIR=$(echo $INPUT_ARGS | cut -d : -f 3)
AREA=$(echo $INPUT_ARGS | cut -d : -f 4)

#### SECTION 4 - Check parameters
# The area must be contained within the sample name and image directory

if [[ "${SAMPLE}" != *"${AREA}"* ]]; then
  echo "Slide area name is not contained in the Sample name. Exiting."
  exit
fi

if [[ "${IMAGE_DIR}" != *"${AREA}"* ]]; then
  echo "Slide area name is not contained in the Image name. Exiting."
  exit
fi

#### SECTION 5    
# Commands to run count
# Ensure that localcores = the number of cores requested in SECTION 1

spaceranger count --id=${SAMPLE}_count \
--fastqs=$FASTQ_DIR \
--sample=$SAMPLE \
--transcriptome=/data/WHRI-GenomeCentre/Genome/cellranger_5.0_references/human/refdata-gex-GRCh38-2020-A/ \
--image=$IMAGE_DIR \
--slide=V10N16-080 \
--area=$AREA \
--localcores=$NSLOTS \
--localmem=8 \
--jobmode=local
