#!/bin/bash
#
#SBATCH --job-name=extract_biosample_matrix
#SBATCH --error=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/error/quantile_normalization.err
#SBATCH --output=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/output/quantile_normalization.out
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=481280

module load anaconda
module load bedops
eval "$(conda shell.bash hook)"
conda activate example_env

# Run Python File
echo "Quantile Normalization Starting"

python /work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/makefiles/quartile_normalization_LW20230123.py

echo "Done"
