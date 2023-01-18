#!/bin/bash
#
#SBATCH --job-name=biosample_selections
#SBATCH --error=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/error/biosamples_selections.err
#SBATCH --output=/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/output/biosamples_selections.out
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=481280

#OOE
echo "Biosample Selection"
<<<<<<< HEAD:region_narration/code/pipeline/ARCHS4_LW/python_files/makefiles/biosample_selection.sh
python /home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/python_files/makefiles/biosamples_selections_LW20230117.py $NEW_MAPPED_FILE $EXPECTED $NUM_STATES $BASE_HDF5_FILE
=======
python /home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/biosamples_selections_LW20230117.py $NEW_MAPPED_FILE $EXPECTED $NUM_STATES $BASE_HDF5_FILE
>>>>>>> b8381f4b579966196fabb88c55dbb3da8663b3a2:region_narration/code/pipeline/ARCHS4_LW/sbatch_files/makefiles/biosample_selection.sh

echo "Done"
