#!/bin/bash
#SBATCH --job-name=DHS_Precision_Recall_CT20211020
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30000
#SBATCH --error=/home/ctrader/Hold/DHS/slurm_output/job.%J.err                          
#SBATCH --output=/home/ctrader/Hold/DHS/slurm_output/job.%J.out
pwd; hostname; date

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/AG/AugmentedGenomics/recommendersystems/lightfm/DHS/

echo "Running DHS Test Precision Recall"
python DHS_precision_recall_CT20211020.py DHS_model_e15_d16_reg0001 DHS_data_e15_d16_reg0001 'test' [1,10,100,1000,10000]
echo "Running DHS Train Precision Recall"
python DHS_precision_recall_CT20211020.py DHS_model_e15_d16_reg0001 DHS_data_e15_d16_reg0001 'train' [1,10,100,1000,10000]
echo "Done running"
