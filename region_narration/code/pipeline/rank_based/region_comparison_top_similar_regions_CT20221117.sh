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

DHS_HDF5_FILE=./regions_25kb_CT20221104.DHSsigmoideuc.hdf5
TF_HDF5_FILE=./regions_25kb_CT20221104.TFsigmoideuc.hdf5
EPILOGOS_HDF5_FILE=./regions_25kb_CT20221104.epilogossigmoideuc.hdf5
COEX_HDF5_FILE=./regions_25kb_CT20221031.coexpression_human_transferred_maxSTRINGsimilarity.hdf5
COOCCUR_HDF5_FILE=./regions_25kb_CT20221031.cooccurenceSTRINGsimilarity.hdf5
DATABASE_HDF5_FILE=./regions_25kb_CT20221031.database_human_transferred_maxSTRINGsimilarity.hdf5
EXPER_HDF5_FILE=./regions_25kb_CT20221031.experiments_human_transferred_maxSTRINGsimilarity.hdf5
FUSION_HDF5_FILE=./regions_25kb_CT20221031.fusionSTRINGsimilarity.hdf5
HOMOLOGY_HDF5_FILE=./regions_25kb_CT20221031.homologySTRINGsimilarity.hdf5
NEIGH_HDF5_FILE=./regions_25kb_CT20221031.neighborhood_human_transferred_maxSTRINGsimilarity.hdf5
TEXTMIN_HDF5_FILE=./regions_25kb_CT20221031.textmining_human_transferred_maxSTRINGsimilarity.hdf5

##########################

#find similar regions
#echo "Bar"
time python /Users/ctrader/Documents/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_top_similar_regions_CT20221117.py $DHS_HDF5_FILE \
$TF_HDF5_FILE $EPILOGOS_HDF5_FILE $COEX_HDF5_FILE $COOCCUR_HDF5_FILE $DATABASE_HDF5_FILE $EXPER_HDF5_FILE $FUSION_HDF5_FILE $HOMOLOGY_HDF5_FILE $NEIGH_HDF5_FILE $TEXTMIN_HDF5_FILE \
$SCALE_GZ_FILENAME $GENOMIC_LOCATION $SEARCH_METHOD
