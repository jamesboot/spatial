#!/bin/bash

# Author: James Boot, Date: 01/07/2024
# Script for running FICTURE on CosMx pre-processed data
# Run format1.sh before this script

#####
#$ -M j.boot@qmul.ac.uk
#$ -m bes
#$ -cwd
#$ -pe smp 4
#$ -l h_vmem=8G  
#$ -j y
#$ -l h_rt=240:00:00
#$ -N fictureRun
#####

# Load Python module
module load python/3.11.6

# Activate virtualenv
source /data/home/hmy961/ficture/bin/activate

# Define input folder
path=/data/home/hmy961/ficture_trial/tma2/outs

# Define scales etc.
mu_scale=1                                          # If your data's coordinates are already in micrometer
key=Count
major_axis=Y                                        # If your data is sorted by the Y-axis

# Create pixel mini batches 
batch_size=500
batch_buff=30
input=${path}/filtered.matrix.TMA2.tsv.gz
batch=${path}/batched.matrix.tsv

ficture make_spatial_minibatch --input ${input} \
--output ${batch} \
--mu_scale ${mu_scale} \
--batch_size ${batch_size} \
--batch_buff ${batch_buff} \
--major_axis ${major_axis}

sort -S 4G -k2,2n -k1,1g ${batch} | gzip -c > ${batch}.gz
rm ${batch}

# Prepare training mini batches
train_width=12                                      # \sqrt{3} x the side length of the hexagon (um)
min_ct_per_unit=50
input=${path}/filtered.matrix.TMA2.tsv.gz
out=${path}/hexagon.d_${train_width}.tsv

ficture make_dge --key ${key} \
--count_header ${key} \
--input ${input} \
--output ${out} \
--hex_width ${train_width} \
--n_move 2 \
--min_ct_per_unit ${min_ct_per_unit} \
--mu_scale ${mu_scale} \
--precision 2 \
--major_axis ${major_axis}

sort -S 4G -k1,1n ${out} | gzip -c > ${out}.gz      # Shuffle hexagons
rm ${out}

# Model training

# Parameters for initializing model
nFactor=12                                          # Number of factors
sliding_step=2
train_nEpoch=3
# train_width=12                                    # Should be the same as used in the above step
model_id=nF${nFactor}.d_${train_width}              # An identifier kept in output file names
min_ct_per_feature=20                               # Ignore genes with total count \< 20
R=10                                                # We use R random initializations and pick one to fit the full model
thread=${NSLOTS}                                    # Number of threads to use

# Model initialization

# Parameters
min_ct_per_unit_fit=20
cmap_name="turbo"

# Output identifiers
model_id=nF${nFactor}.d_${train_width}
output_id=${model_id}
output_path=${path}/analysis/${model_id}
figure_path=${output_path}/figure
if [ ! -d "${figure_path}/sub" ]; then
    mkdir -p ${figure_path}/sub
fi

# Input files
hexagon=${path}/hexagon.d_${train_width}.tsv.gz
pixel=${path}/filtered.matrix.TMA2.tsv.gz
feature=${path}/feature.clean.TMA2.tsv.gz

# Output
output=${output_path}/${output_id}
model=${output}.model.p

# Fit model
ficture fit_model --input ${hexagon} \
--output ${output} \
--feature ${feature} \
--nFactor ${nFactor} \
--epoch ${train_nEpoch} \
--epoch_id_length 2 \
--unit_attr X Y \
--key ${key} \
--min_ct_per_feature ${min_ct_per_feature} \
--test_split 0.5 \
--R ${R} \
--thread ${thread}

# Choose color
input=${output_path}/${output_id}.fit_result.tsv.gz
output=${figure_path}/${output_id}
cmap=${figure_path}/${output_id}.rgb.tsv

ficture choose_color --input ${input} \
--output ${output} \
--cmap_name ${cmap_name}

# Coarse plot for inspection
cmap=${figure_path}/${output_id}.rgb.tsv
input=${output_path}/${output_id}.fit_result.tsv.gz
output=${figure_path}/${output_id}.coarse
fillr=$((train_width/2+1))

ficture plot_base --input ${input} \
--output ${output} \
--fill_range ${fillr} \
--color_table ${cmap} \
--plot_um_per_pixel 1 \
--plot_discretized

# Pixel level decoding

# Parameters
fit_width=12                                                            # Often equal or smaller than train_width (um)
anchor_res=4                                                            # Distance between adjacent anchor points (um)
fit_nmove=$((fit_width/anchor_res))
anchor_info=prj_${fit_width}.r_${anchor_res}
radius=$(($anchor_res+1))
anchor_info=prj_${fit_width}.r_${anchor_res}                            # An identifier
coor=${path}/coordinate_minmax.tsv
cmap=${figure_path}/${output_id}.rgb.tsv

# Transform
output=${output_path}/${output_id}.${anchor_info}
ficture transform --input ${pixel} \
--output_pref ${output} \
--model ${model} \
--key ${key} \
--major_axis ${major_axis} \
--hex_width ${fit_width} \
--n_move ${fit_nmove} \
--min_ct_per_unit ${min_ct_per_unit_fit} \
--mu_scale ${mu_scale} \
--thread ${thread} \
--precision 2

# Pixel level decoding & visualization
prefix=${output_id}.decode.${anchor_info}_${radius}
input=${path}/batched.matrix.tsv.gz
anchor=${output_path}/${output_id}.${anchor_info}.fit_result.tsv.gz
output=${output_path}/${prefix}
topk=5                                                                  # Output only a few top factors per pixel

ficture slda_decode --input ${input} \
--output ${output} \
--model ${model} \
--anchor ${anchor} \
--anchor_in_um \
--neighbor_radius ${radius} \
--mu_scale ${mu_scale} \
--key ${key} \
--precision 0.1 \
--lite_topk_output_pixel ${topk} \
--lite_topk_output_anchor ${topk} \
--thread ${thread}

# Optional post-processing
input=${output_path}/${prefix}.pixel.tsv.gz                             # j, X, Y, K1, ..., KJ, P1, ..., PJ, J=topk
output=${output_path}/${prefix}.pixel.sorted.tsv.gz

K=$( echo $model_id | sed 's/nF\([0-9]\{1,\}\)\..*/\1/' )
while IFS=$'\t' read -r r_key r_val; do
    export "${r_key}"="${r_val}"
done < ${coor}
echo -e "${xmin}, ${xmax}; ${ymin}, ${ymax}"

offsetx=${xmin}
offsety=${ymin}
rangex=$( echo "(${xmax} - ${xmin} + 0.5)/1+1" | bc )
rangey=$( echo "(${ymax} - ${ymin} + 0.5)/1+1" | bc )
bsize=2000
scale=100
header="##K=${K};TOPK=3\n##BLOCK_SIZE=${bsize};BLOCK_AXIS=X;INDEX_AXIS=Y\n##OFFSET_X=${offsetx};OFFSET_Y=${offsety};SIZE_X=${rangex};SIZE_Y=${rangey};SCALE=${scale}\n#BLOCK\tX\tY\tK1\tK2\tK3\tP1\tP2\tP3"

(echo -e "${header}" && zcat ${input} | tail -n +2 | perl -slane '$F[0]=int(($F[1]-$offx)/$bsize) * $bsize; $F[1]=int(($F[1]-$offx)*$scale); $F[1]=($F[1]>=0)?$F[1]:0; $F[2]=int(($F[2]-$offy)*$scale); $F[2]=($F[2]>=0)?$F[2]:0; print join("\t", @F);' -- -bsize=${bsize} -scale=${scale} -offx=${offsetx} -offy=${offsety} | sort -S 4G -k1,1g -k3,3g ) | bgzip -c > ${output}

tabix -f -s1 -b3 -e3 ${output}
rm ${input}

# DE Analysis
max_pval_output=1e-3
min_fold_output=1.5
input=${output_path}/${prefix}.posterior.count.tsv.gz
output=${output_path}/${prefix}.bulk_chisq.tsv

ficture de_bulk --input ${input} \
--output ${output} \
--min_ct_per_feature ${min_ct_per_feature} \
--max_pval_output ${max_pval_output} \
--min_fold_output ${min_fold_output} \
--thread ${thread}

# Report (color table and top DE genes)
cmap=${output_path}/figure/${output_id}.rgb.tsv
output=${output_path}/${prefix}.factor.info.html

ficture factor_report --path ${output_path} \
--pref ${prefix} \
--color_table ${cmap}

# Make pixel level figures
cmap=${output_path}/figure/${output_id}.rgb.tsv
input=${output_path}/${prefix}.pixel.sorted.tsv.gz
output=${figure_path}/${prefix}.pixel.png

ficture plot_pixel_full --input ${input} \
--color_table ${cmap} \
--output ${output} \
--plot_um_per_pixel 0.5 \
--full

# Make single factor heatmaps
output=${figure_path}/sub/${prefix}.pixel

ficture plot_pixel_single --input ${input} \
--output ${output} \
--plot_um_per_pixel 0.5 \
--full \
--all