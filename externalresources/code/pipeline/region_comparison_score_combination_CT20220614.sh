#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/score_comb_job.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/score_comb_job.out
pwd; hostname; date

module load anaconda
#module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data

#####Variables####
SORTED_SCALE_FILENAME=Altius_25kb_hg38_maxMean_CT20220613.sorted.c1t3.bed
SCALE_GZ_FILENAME=Altius_25kb_hg38_maxMean_CT20220613.sorted.c1t3.rownum4.bed.gz

NUM_CPUS=1

DHS_COSINE_HDF5_FILE=./regions_25kb_CT20220613.DHScosine.hdf5
TF_COSINE_HDF5_FILE=./regions_25kb_CT20220613.TFcosine.hdf5
EPILOGOS_COSINE_HDF5_FILE=./regions_25kb_CT20220613.epilogoscosine.hdf5
##################
GENOMIC_LOCATION=$1
##################

echo $GENOMIC_LOCATION

#find similar regions
echo "Finding Similar Regions"
time python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_combine_scores_CT20220614.py $SCALE_GZ_FILENAME $GENOMIC_LOCATION $DHS_COSINE_HDF5_FILE $TF_COSINE_HDF5_FILE $EPILOGOS_COSINE_HDF5_FILE $SORTED_SCALE_FILENAME

echo "Done"
