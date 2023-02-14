#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/signal_preprocessing.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/signal_preprocessing.out
pwd; hostname; date

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data/rank_based

#####Variables####

#DHS
MAPPED_FILE=regions_25kb_CT20230118.DHSsignalmapped.bed
Z_HDF5_FILE=regions_25kb_CT20230202.DHSsignal.factorblindzscores.hdf5

#TF Motifs
#MAPPED_FILE=regions_25kb_CT20230118.TFsignalmapped.bed
#Z_HDF5_FILE=regions_25kb_CT20230118.TFsignalzscores.hdf5
##################

#factor blind z score of signal data
echo "Factor blind Z Score"
python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_signal_factor_blind_zscore_CT20230202.py  $MAPPED_FILE $Z_HDF5_FILE

echo "Done"
