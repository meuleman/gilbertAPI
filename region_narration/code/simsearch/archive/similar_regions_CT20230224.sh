#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate simsearch_env

#data dir
DATADIR=/Users/ctrader/Documents/local_data/region_similarity/data/rank_based/
SIMSEARCHDIR=/Users/ctrader/Documents/AG/NarrativeGenomics/region_narration/code/simsearch/
##########################

GENOMIC_LOCATION=$1
TARGET_FACTORS=$2
AUTO_FACTOR_SELECTION=$3

##########################
#data
SCALE_GZ_FILENAME=${DATADIR}hg38.25kb.maxMean.CT20230118.sorted.rownum4.bed.gz
REGIONS_PATH=${DATADIR}hg38.25kb.maxMean.CT20230118.sorted.bed

DHS_HDF5_FILE=${DATADIR}regions_25kb_CT20230214.DHS.factornormNMF.sum_enrich_corank_percentile.hdf5
EPILOGOS_HDF5_FILE=${DATADIR}regions_25kb_CT20230214.epilogos.statesum_rank_percentile.hdf5
##########################


#find similar regions
time python ${SIMSEARCHDIR}similar_regions_CT20230224.py $DHS_HDF5_FILE \
$EPILOGOS_HDF5_FILE $SCALE_GZ_FILENAME $REGIONS_PATH $GENOMIC_LOCATION $TARGET_FACTORS \
$AUTO_FACTOR_SELECTION
