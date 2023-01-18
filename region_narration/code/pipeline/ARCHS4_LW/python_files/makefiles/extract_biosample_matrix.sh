#!/bin/bash
#
#SBATCH --job-name=extract_biosample_matrix
#SBATCH --error=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/python_files/error/extract_biosample_matrix.err
#SBATCH --output=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/python_files/output/extract_biosample_matrix.out
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=481280

# Run Python File
echo "Extracting Biosample Filtered Matrix"
python /home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/python_files/makefiles/extract_biosample_matrix.py

echo "Done"
