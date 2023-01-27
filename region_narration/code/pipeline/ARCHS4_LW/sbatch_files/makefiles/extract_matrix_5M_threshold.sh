#!/bin/bash
#
#SBATCH --job-name=extract_matrix_5M
#SBATCH --error=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/error/extract_matrix_5M.err
#SBATCH --output=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/output/extract_matrix_5M.out
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=481280

module load anaconda
module load bedops
eval "$(conda shell.bash hook)"
conda activate example_env

# Run Python File
echo "Extracting 5M Threshold Matrix"

python /home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/makefiles/extract_matrix_5M_threshold_LW20230125.py

echo "Done"
