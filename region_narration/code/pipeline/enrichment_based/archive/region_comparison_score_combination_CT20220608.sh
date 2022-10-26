#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/job.%J.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/job.%J.out

#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/score_comb_job.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/score_comb_job.out
pwd; hostname; date

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data

#####Variables####
SORTED_SCALE_FILENAME=Altius_100k_epilogos_hg19_WM20210518.sorted.bed
REGION_FILENAME=./region_of_interest_CT20220606.txt
NUM_CPUS=1

DHS_COSINE_NPZ_FILE=./regions_100k_CT20220602.DHScosine.npz
TF_COSINE_NPZ_FILE=./regions_100k_CT20220602.TFcosine.npz
EPILOGOS_COSINE_NPZ_FILE=./regions_100k_CT20220602.epilogoscosine.npz
##################
GENOMIC_LOCATION=$1
##################

echo $GENOMIC_LOCATION

#find region
echo "Finding Region of Interest"
python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_find_region_ind_CT20220606.py $SORTED_SCALE_FILENAME $GENOMIC_LOCATION $REGION_FILENAME

#find similar regions
echo "Finding Similar Regions"
python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_combine_scores_CT20220609.py $REGION_FILENAME $DHS_COSINE_NPZ_FILE $TF_COSINE_NPZ_FILE $EPILOGOS_COSINE_NPZ_FILE

echo "Done"
