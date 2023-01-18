#!/bin/bash
#
#SBATCH --job-name=biosample_selections
#SBATCH --error=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/error/biosamples_selections.err
#SBATCH --output=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/output/biosamples_selections.out
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=481280

# Run Python File
echo "Biosample Selection"

python /home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/biosamples_selections_LW20230117.py $NEW_MAPPED_FILE $EXPECTED $NUM_STATES $BASE_HDF5_FILE

echo "Done"
