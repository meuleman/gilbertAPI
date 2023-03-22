#!/bin/bash

genomic_location=${1}
included_factors=${2}
get_roiURL=${3}
alternate_data_dir=${4}
alternate_simsearch_dir=${5}

if [ ! -n "${genomic_location}" ]
then
    echo "Error: Genomic location was not provided"
    exit -1
fi

if [ ! -n "${included_factors}" ]
then
    echo "Error: Factors were not provided"
    exit -1
fi

if [ ! -n "${get_roiURL}" ]
then
    echo "Error: option to get roiURL not selected"
    exit -1
fi

##########################
# Allow omission to replace with useful defaults
##########################

if [ ! -n "${alternate_data_dir}" ]
then
    data_dir=/Users/ctrader/Documents/AG/NarrativeGenomics/region_narration/data
else
    data_dir="${alternate_data_dir}"
fi

if [ ! -n "${alternate_simsearch_dir}" ]
then
    simsearch_dir=/Users/ctrader/Documents/AG/NarrativeGenomics/region_narration/code/simsearch
else
    simsearch_dir="${alternate_simsearch_dir}"
fi

##########################
# Altius data dirs
# DATADIR=/Users/ctrader/Documents/local_data/region_similarity/data/rank_based
# SIMSEARCHDIR=/Users/ctrader/Documents/AG/NarrativeGenomics/region_narration/code/simsearch
altius_data_dir=/home/ctrader/AG/NarrativeGenomics/region_narration/data
altius_simsearch_dir=/home/ctrader/AG/NarrativeGenomics/region_narration/code/simsearch
##########################

##########################
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
if [ -d "${data_dir}" ]
then
    true_data_dir="${data_dir}"
    true_code_dir="${simsearch_dir}"
fi

scale_gz_filename="${true_data_dir}/hg38.25kb.maxMean.CT20230118.sorted.rownum4.bed.gz"
regions_filename="${true_data_dir}/hg38.25kb.maxMean.CT20230118.sorted.bed.npy"
per_filename="${true_data_dir}/regions_25kb_CT20230321.DHS.factornormNMF.sum_enrich_corank.epilogos.statesum_rank.percentiles.npy"
num_similar_regions=20
##########################

##########################
# find similar regions
#
# using cProfile to profile function runtimes
#
#  $ pip install cprofilev (if needed)
#  $ cprofilev -f ${region_narration_express_dir}/similar_regions.cProfile.out
#
# open http://127.0.0.1:4000/
#
##########################
time python -m cProfile -o similar_regions.cProfile.out "${true_code_dir}/similar_regions.py" \
    "${per_filename}" \
    "${scale_gz_filename}" \
    "${regions_filename}" \
    "${genomic_location}" \
    "${included_factors}" \
    "${get_roiURL}" \
    "${num_similar_regions}"
##########################
