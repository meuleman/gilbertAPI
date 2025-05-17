#!/bin/bash
#SBATCH --job-name=model_train_array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=18000
#SBATCH --error=/home/ctrader/Hold/GWAS/slurm_output/job.%J.err
#SBATCH --output=/home/ctrader/Hold/GWAS/slurm_output/job.%J.out
#SBATCH --array=[5,10,15,20,25,40,50,75,100]
pwd; hostname; date

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/AG/AugmentedGenomics/recommendersystems/lightfm/GWAS/

echo "Running job array"
python train_GWAS_model_afo_alpha_epochs_dimensions_CT20211130.py 50 $SLURM_ARRAY_TASK_ID 20
echo "Done running"
