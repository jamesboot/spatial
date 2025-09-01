# This script is for renaming fastq files if they are missing L001 from file name
# Needs to be run before GeoMx NGS Pipeline or will result in error

#!/bin/bash

##### Set queue options
#$ -M j.boot@qmul.ac.uk					# Email address
#$ -m bes											# Send email
#$ -cwd
#$ -pe smp 1										# Cores
#$ -l h_vmem=1G								# Memory
#$ -l h_rt=1:00:00							# Running time
#$ -N Rename_FastQ					# Rename this
#$ -j y
#####

# For R1s
for FILE in *_R1_001.fastq.gz; do
	mv ${FILE} $(echo ${FILE} | sed 's/_R1/_L001_R1/g')
done

# For R2s
for FILE in *_R2_001.fastq.gz; do
	mv ${FILE} $(echo ${FILE} | sed 's/_R2/_L001_R2/g')
done

