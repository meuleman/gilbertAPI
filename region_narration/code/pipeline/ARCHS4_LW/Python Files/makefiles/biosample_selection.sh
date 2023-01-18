#!/bin/bash
#
#SBATCH --job-name=biosample_selections
#SBATCH --error=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/Python Files/error/biosamples_selections.err
#SBATCH --output=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/Python Files/output/biosamples_selections.out
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100000

#OOE
echo "Biosample Selection"
python /home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/Python Files/biosamples_selections_LW20230117.py $NEW_MAPPED_FILE $EXPECTED $NUM_STATES $BASE_HDF5_FILE

echo "Done"
