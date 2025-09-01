# This script is for running GeoMx NGS Pipeline to pre-process data and generate DCC files

#!/bin/bash

##### Set queue options
#$ -M j.boot@qmul.ac.uk					# Email address
#$ -m bes											# Send email
#$ -cwd
#$ -pe smp 4										# Cores
#$ -l h_vmem=32G								# Memory
#$ -l h_rt=240:00:00							# Running time
#$ -N GeoMxNGSPipeline					# Rename this
#$ -j y
#####

CONFIG=/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/NanostringGeoMX/Projects/Martin-Gutierrez-Lucia/GC-LM-10326/Analysis/GeoMx_NGS_Pipeline_Outs/GC-LM-10326_20230606T1621_GNP_config.ini
OUT=/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/NanostringGeoMX/Projects/Martin-Gutierrez-Lucia/GC-LM-10326/Analysis/GeoMx_NGS_Pipeline_Outs
IN=/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/NanostringGeoMX/Projects/Martin-Gutierrez-Lucia/GC-LM-10326/Data/230615_VH01330_50_AAC3FCGHV/Analysis/3/Data/Intensities/BaseCalls

/data/WHRI-GenomeCentre/software/GeoMxNGSPipeline/geomxngspipeline --in=$IN \
--out=$OUT \
--ini=$CONFIG
