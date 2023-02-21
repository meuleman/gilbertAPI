#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate lightfm

#data dir
cd /home/ctrader/region_similarity/data/rank_based/

##########################

GENOMIC_LOCATION=$1
TARGET_FACTORS=$2

##########################
#data
SCALE_GZ_FILENAME=hg38.25kb.maxMean.CT20230118.sorted.rownum4.bed.gz
REGIONS_PATH=hg38.25kb.maxMean.CT20230118.sorted.bed

DHS_HDF5_FILE=./regions_25kb_CT20230214.DHS.factornormNMF.sum_enrich_corank_percentile.hdf5
EPILOGOS_HDF5_FILE=./regions_25kb_CT20230214.epilogos.statesum_rank_percentile.hdf5
##########################



#TARGET_FACTORS=[3,30]


#find similar regions
time python /home/ctrader/AG/NarrativeGenomics/region_narration/code/simsearch/similar_regions_roiurl_CT20230215.py $DHS_HDF5_FILE \
$EPILOGOS_HDF5_FILE $SCALE_GZ_FILENAME $REGIONS_PATH $GENOMIC_LOCATION $TARGET_FACTORS

