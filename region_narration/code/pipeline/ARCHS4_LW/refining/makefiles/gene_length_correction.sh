#!/bin/bash
#
#SBATCH --job-name=extract_biosample_matrix
#SBATCH --error=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/error/gene_length_correction.err
#SBATCH --output=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/output/gene_length_correction.out
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=481280

module load anaconda
module load bedops
eval "$(conda shell.bash hook)"
conda activate example_env

# Run Python File
echo "Computing Gene Length Correction Matrix"

python /work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/refining/gene_length_correction/gene_length_correction.py

echo "Done"
