#!/bin/bash
#
#SBATCH --job-name=extract_biosample_matrix
#SBATCH --error=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/python_files/error/biosamples_selections.err
#SBATCH --output=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/python_files/output/biosamples_selections.out
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100000

#OOE
echo "Extracting Biosample Filtered Matrix"
python /home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/python_files/makefiles/extract_biosample_matrix.py

echo "Done"
