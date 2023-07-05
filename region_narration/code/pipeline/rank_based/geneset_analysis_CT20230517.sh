#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/geneset.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/geneset.out
pwd; hostname; date
today=$(date "+%Y%m%d")

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

python geneset_analysis_CT20230517.py
