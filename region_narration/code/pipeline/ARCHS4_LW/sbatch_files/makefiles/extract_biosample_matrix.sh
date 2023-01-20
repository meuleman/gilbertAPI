#!/bin/bash
#
#SBATCH --job-name=extract_biosample_matrix
#SBATCH --error=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/error/biosamples_selections.err
#SBATCH --output=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/output/biosamples_selections.out
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=481280

module load anaconda
conda activate example_env

# Run Python File
echo "Extracting Biosample Filtered Matrix"

python /home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/makefiles/extract_biosample_matrix.py

echo "Done"
