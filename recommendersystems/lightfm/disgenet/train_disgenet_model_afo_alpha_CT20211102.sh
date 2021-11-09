#!/bin/bash
#SBATCH --job-name=model_train_array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=18000
#SBATCH --error=/home/ctrader/Hold/Dis/slurm_output/job.%J.err
#SBATCH --output=/home/ctrader/Hold/Dis/slurm_output/job.%J.out
#SBATCH --array=[0,1,5,10,50,100,500,1000]
pwd; hostname; date

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/AG/AugmentedGenomics/recommendersystems/lightfm/disgenet/

echo "Running job array"
python train_disgenet_model_afo_alpha_CT20211102.py $SLURM_ARRAY_TASK_ID
echo "Done running"
