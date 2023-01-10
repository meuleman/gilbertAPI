#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate region_rec_env

cd /Users/ctrader/Documents/local_data/region_similarity/data

##########################

GENOMIC_LOCATION=$1
SEARCH_METHOD=$2

####Variable to change####

#FIGURE_OUT_DIR=/Users/ctrader/Documents/local_data/region_similarity/data/local_figures
#SENTENCE_OUT=/Users/ctrader/Documents/local_data/region_similarity/data/sentence_summaries/regions_25kb_${GENOMIC_LOCATION/:/_}_sentence_summaries_CT$(date +%Y%m%d).txt
#UPSET_LIMITING_PER=95

##########################

SCALE_GZ_FILENAME=Altius_25kb_hg38_maxMean_CT20220815.sorted.c1t3.rownum4.bed.gz
REGIONS_PATH=Altius_25kb_hg38_maxMean_CT20220815.sorted.c1t3.bed

DHS_HDF5_FILE=./regions_25kb_CT20221104.DHSsigmoideuc.hdf5
TF_HDF5_FILE=./regions_25kb_CT20221104.TFsigmoideuc.hdf5
EPILOGOS_HDF5_FILE=./regions_25kb_CT20221104.epilogossigmoideuc.hdf5

DHS_MAPPED_FILE=./regions_25kb_CT20221104.DHSmapped.bed
TF_MAPPED_FILE=./regions_25kb_CT20221104.TFmapped.bed
EPILOGOS_MAPPED_FILE=./regions_25kb_CT20221104.epilogosmapped.bed
STRING_RANKS_HDF5_FILE=./regions_25kb_CT20221115.top_gene.combined_scores_above_400.counts_by_metric.hdf5
##########################

#find similar regions
#echo "Bar"
time python /Users/ctrader/Documents/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_top_similar_regions_minimum_CT20230110.py $DHS_HDF5_FILE \
$TF_HDF5_FILE $EPILOGOS_HDF5_FILE $DHS_MAPPED_FILE $TF_MAPPED_FILE $EPILOGOS_MAPPED_FILE $STRING_RANKS_HDF5_FILE $SCALE_GZ_FILENAME $REGIONS_PATH $GENOMIC_LOCATION $SEARCH_METHOD

