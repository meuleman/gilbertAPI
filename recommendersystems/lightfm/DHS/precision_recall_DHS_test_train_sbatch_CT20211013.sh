#!/bin/bash
#SBATCH --job-name=DHS_Precision_Recall_CT20211013
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=20000
#SBATCH --error=/home/ctrader/AG/Hold/job.%J.err
#SBATCH --output=/home/ctrader/AG/Hold/job.%J.out
pwd; hostname; date

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/AG/AugmentedGenomics/recommendersystems/lightfm/DHS/

echo "Running DHS_Test_Precision_Recall"
python precision_recall_function_DHS_test_CT20211013.py
echo "Running DHS_Train_Precision_Recall"
python precision_recall_function_DHS_train_CT20211013.py
echo "Done running"
