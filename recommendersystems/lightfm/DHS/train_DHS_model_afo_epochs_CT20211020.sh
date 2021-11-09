#!/bin/bash
#SBATCH --job-name=model_train_array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=18000
#SBATCH --error=/home/ctrader/Hold/DHS/slurm_output/job.%J.err
#SBATCH --output=/home/ctrader/Hold/DHS/slurm_output/job.%J.out
#SBATCH --array=[0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009]
pwd; hostname; date

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/AG/AugmentedGenomics/recommendersystems/lightfm/DHS/

echo "Running job array"
python train_DHS_model_afo_lr_CT20211026.py $SLURM_ARRAY_TASK_ID
echo "Done running"
