#!/bin/bash

GENOMIC_LOCATION=${1}
INCLUDED_FACTORS=${2}
GET_ROIURL=${3}

if [ ! -n "${GENOMIC_LOCATION}" ]
then
    echo "Error: Genomic location was not provided"
    exit -1
fi

if [ ! -n "${INCLUDED_FACTORS}" ]
then
    echo "Error: Factors were not provided"
    exit -1
fi
if [ ! -n "${GET_ROIURL}" ]
then
    echo "Error: option to get roiURL not selected"
    exit -1
fi

##########################
# data dir
DATADIR=/Users/ctrader/Documents/local_data/region_similarity/data/rank_based
altius_data_dir=/home/ctrader/AG/NarrativeGenomics/region_narration/data
SIMSEARCHDIR=/Users/ctrader/Documents/AG/NarrativeGenomics/region_narration/code/simsearch
altius_simsearch_dir=/home/ctrader/AG/NarrativeGenomics/region_narration/code/simsearch
##########################

#########################
# test if conda is available
environment_name=simsearch_env
if command -v conda &> /dev/null
then
    eval "$(conda shell.bash hook)"
    conda activate ${environment_name}
else
    requirements_fn=${altius_simsearch_dir}/requirements.txt
    module add anaconda
    eval "$(conda shell.bash hook)"
    conda env list
    conda config --append channels conda-forge
    conda config --append channels bioconda
    conda create --name ${environment_name} python=3.8 -y
    conda activate ${environment_name}
    pip install -r ${requirements_fn}
fi

##########################
# data
true_data_dir=${altius_data_dir}
true_code_dir=${altius_simsearch_dir}
if [ -d ${DATADIR} ]
then
    true_data_dir=${DATADIR}
    true_code_dir=${SIMSEARCHDIR}
fi
SCALE_GZ_FILENAME=${true_data_dir}/hg38.25kb.maxMean.CT20230118.sorted.rownum4.bed.gz
REGIONS_PATH=${true_data_dir}/hg38.25kb.maxMean.CT20230118.sorted.bed
DHS_HDF5_FILE=${true_data_dir}/regions_25kb_CT20230214.DHS.factornormNMF.sum_enrich_corank_percentile.hdf5
EPILOGOS_HDF5_FILE=${true_data_dir}/regions_25kb_CT20230214.epilogos.statesum_rank_percentile.hdf5
##########################

##########################
# find similar regions
time python ${true_code_dir}/similar_regions.py \
    $DHS_HDF5_FILE \
    $EPILOGOS_HDF5_FILE \
    $SCALE_GZ_FILENAME \
    $REGIONS_PATH \
    $GENOMIC_LOCATION \
    $INCLUDED_FACTORS \
    $GET_ROIURL
##########################
