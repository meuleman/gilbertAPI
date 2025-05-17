#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/genesetXgene_clustering.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/genesetXgene_clustering.out
pwd; hostname; date
today=$(date "+%Y%m%d")

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

python genesetXgene_clustering_CT20230606.py
