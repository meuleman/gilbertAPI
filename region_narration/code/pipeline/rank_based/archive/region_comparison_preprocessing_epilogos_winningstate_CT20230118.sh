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
WINNING_STATES=epilogos.winningstateperbin.officialorder.CT20230118.sorted.bed
SORTED_SCALE_FILENAME=hg38.25kb.maxMean.CT20230118.sorted.bed
EXPECTED=epilogos.winningstateperbin.genomewide.counts.expectedratio.officialorder.CT20230118.tsv
NUM_STATES=18

#Variables to change
NEW_MAPPED_FILE=regions_25kb_CT20230118.epilogos.maxstateperbin.officialstateorder.allmapped.bed
BASE_HDF5_FILE=regions_25kb_CT20230118.epilogos.maxstateperbin.officialstateorder.ooe.hdf5
CONV_HDF5_FILE=regions_25kb_CT20230118.epilogos.maxstateperbin.officialstateorder.rank.sigmoid.conv_sqrt_sigmoid_sim.max_priors.hdf5
##################

#mapping winning states to sorted scale file
bedmap --echo --echo-map-score $SORTED_SCALE_FILENAME $WINNING_STATES > $NEW_MAPPED_FILE

#OOE
echo "OOE"
python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_epilogos_winningstate_ooe_CT20230112.py $NEW_MAPPED_FILE $EXPECTED $NUM_STATES $BASE_HDF5_FILE

#Sigmoid Euclidean Similarity and Priors
echo "Calculating Similarity and Priors"
python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_rank_sigmoid_convsim_priors_CT20230117.py $BASE_HDF5_FILE $CONV_HDF5_FILE

echo "Done"
