#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/score_comb_job.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/score_comb_job.out
pwd; hostname; date

module load anaconda
#module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data

#####Variables####
SORTED_SCALE_FILENAME=Altius_100k_epilogos_hg19_WM20210518.sorted.c1t3.bed
#REGION_FILENAME=./region_of_interest_CT20220606.txt
NUM_CPUS=1

DHS_COSINE_HDF5_FILE=./regions_100k_CT20220610.DHScosine.hdf5
TF_COSINE_HDF5_FILE=./regions_100k_CT20220610.TFcosine.hdf5
EPILOGOS_COSINE_HDF5_FILE=./regions_100k_CT20220610.epilogoscosine.hdf5
##################
GENOMIC_LOCATION=$1
##################

echo $GENOMIC_LOCATION

#find region
echo "Finding Region of Interest"
#time python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_find_region_ind_CT20220606.py $SORTED_SCALE_FILENAME $GENOMIC_LOCATION $REGION_FILENAME
time REGION_IND=$(tabix Altius_100k_epilogos_hg19_WM20210518.sorted.c1t3.1row4.bed.gz $GENOMIC_LOCATION | cut -f 4)
REGION_IND=$((REGION_IND-1))

echo $REGION_IND

#find similar regions
echo "Finding Similar Regions"
time python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_combine_scores_CT20220613.py $REGION_IND $DHS_COSINE_HDF5_FILE $TF_COSINE_HDF5_FILE $EPILOGOS_COSINE_HDF5_FILE $SORTED_SCALE_FILENAME

echo "Done"
