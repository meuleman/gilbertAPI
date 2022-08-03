#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/job.%J.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/job.%J.out

#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/job.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/job.out
pwd; hostname; date

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data

#####Variables####
PROMOTER_MAPPED=./regions_25kb_CT20220726.genepromoterbpmapped.bed

NEW_MAPPED_FILE=regions_25kb_CT20220708.epilogosmapped.bed
OOE_HDF5_FILE=regions_25kb_CT20220727.epilogospromoteronlyoee.hdf5
##################

#Calculating OOE
echo "Calculating OOE"
python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_ooe_epilogos_promoteronly_CT20220727.py $NEW_MAPPED_FILE $PROMOTER_MAPPED $OOE_HDF5_FILE

echo "Done"
