#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/preprocessing_epilogos.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/preprocessing_epilogos.out
pwd; hostname; date

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data/rank_based

#####Variables####
WINNING_STATES=../epilogos_scores_new.winningstate.officialorder.sorted.CT20230103.bed
SORTED_SCALE_FILENAME=../Altius_25kb_hg38_maxMean_CT20220815.sorted.c1t3.bed
EXPECTED=epilogos.winningstateperbin.genomewide.counts.expectedratio.officialorder.CT20230103.tsv
NUM_STATES=18

#Variables to change
NEW_MAPPED_FILE=regions_25kb_CT20230112.epilogos.maxstateperbin.officialstateorder.allmapped.bed
BASE_HDF5_FILE=regions_25kb_CT20230112.epilogos.maxstateperbin.officialstateorder.ooe.hdf5
CONV_HDF5_FILE=regions_25kb_CT20230112.epilogos.maxstateperbin.officialstateorder.rank.sigmoid.conv.sqrt.sigmoid.hdf5
##################

#mapping winning states to sorted scale file
bedmap --echo --echo-map-score $SORTED_SCALE_FILENAME $WINNING_STATES > $NEW_MAPPED_FILE

#OOE
echo "OOE"
python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_epilogos_winningstate_ooe_CT20230112.py $NEW_MAPPED_FILE $EXPECTED $NUM_STATES $BASE_HDF5_FILE

#Sigmoid Euclidean Similarity
echo "Calculating Similarity"
python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_rank_sigmoid_conv_sim_CT20230112.py $BASE_HDF5_FILE $CONV_HDF5_FILE

echo "Done"
