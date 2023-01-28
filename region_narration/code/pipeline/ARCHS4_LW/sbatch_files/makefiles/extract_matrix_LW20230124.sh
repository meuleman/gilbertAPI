#!/bin/bash
#
#SBATCH --job-name=extract_matrix
#SBATCH --error=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/error/extract_matrix.err
#SBATCH --output=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/output/extract_matrix.out
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=481280

module load anaconda
module load bedops
eval "$(conda shell.bash hook)"
conda activate example_env

# Run Python File
echo "Extracting Filtered Matrix"

python /home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/makefiles/extract_matrix_LW20230124.py

echo "Done"
