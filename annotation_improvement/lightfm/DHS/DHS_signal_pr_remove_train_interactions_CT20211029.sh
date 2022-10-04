#!/bin/bash
#SBATCH --job-name=DHS_PR_CT20211029
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30000
#SBATCH --error=/home/ctrader/Hold/DHS/pr_out/job.%J.err
#SBATCH --output=/home/ctrader/Hold/DHS/pr_out/job.%J.out
pwd; hostname; date

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/AG/AugmentedGenomics/recommendersystems/lightfm/DHS/

echo "Running DHS Test Precision Recall"
python DHS_signal_pr_remove_train_interactions_CT20211029.py DHS_signal_model_e15_d16_reg0.0001_lr01_rs3 DHS_signal_data_e15_d16_reg0.0001_lr01_rs3 'test'
echo "Running DHS Train Precision Recall"
python DHS_signal_pr_remove_train_interactions_CT20211029.py DHS_signal_model_e15_d16_reg0.0001_lr01_rs3 DHS_signal_data_e15_d16_reg0.0001_lr01_rs3 'train'
echo "Done running"
