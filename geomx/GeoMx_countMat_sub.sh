#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_rt=24:0:0
#$ -l h_vmem=8G

module load R
Rscript GeoMx_countMat.R