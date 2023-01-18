#!/bin/bash
#
#SBATCH --job-name=extract_biosample_matrix
<<<<<<< HEAD:region_narration/code/pipeline/ARCHS4_LW/python_files/makefiles/extract_biosample_matrix.sh
#SBATCH --error=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/python_files/error/extract_biosample_matrix.err
#SBATCH --output=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/python_files/output/extract_biosample_matrix.out
=======
#SBATCH --error=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/error/biosamples_selections.err
#SBATCH --output=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/output/biosamples_selections.out
>>>>>>> b8381f4b579966196fabb88c55dbb3da8663b3a2:region_narration/code/pipeline/ARCHS4_LW/sbatch_files/makefiles/extract_biosample_matrix.sh
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=481280

# Run Python File
echo "Extracting Biosample Filtered Matrix"
python /home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/makefiles/extract_biosample_matrix.py

echo "Done"
