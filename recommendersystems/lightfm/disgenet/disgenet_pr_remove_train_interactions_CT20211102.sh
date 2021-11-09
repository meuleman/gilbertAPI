#!/bin/bash
#SBATCH --job-name=disgenet_PR_CT20211102
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30000
#SBATCH --error=/home/ctrader/Hold/Dis/pr_out/job.%J.err
#SBATCH --output=/home/ctrader/Hold/Dis/pr_out/job.%J.out
pwd; hostname; date

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/AG/AugmentedGenomics/recommendersystems/lightfm/disgenet/

echo "Running disgenet Test Precision Recall"
python disgenet_pr_remove_train_interactions_CT20211102.py disgenet_model_e100_d10_reg5e-05_lr05 disgenet_data_e100_d10_reg5e-05_lr05 'test'
echo "Running disgenet Train Precision Recall"
python disgenet_pr_remove_train_interactions_CT20211102.py disgenet_model_e100_d10_reg5e-05_lr05 disgenet_data_e100_d10_reg5e-05_lr05 'train'
echo "Done running"
