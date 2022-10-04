#!/bin/bash
#SBATCH --job-name=GWAS_PR_CT20211130
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90000
#SBATCH --error=/home/ctrader/Hold/GWAS/pr_out/job.%J.err
#SBATCH --output=/home/ctrader/Hold/GWAS/pr_out/job.%J.out
pwd; hostname; date

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/AG/AugmentedGenomics/recommendersystems/lightfm/GWAS/

echo "Running DHS Tune Precision Recall"
python GWAS_pr_remove_train_interactions_CT20211130.py GWAS_model_e30_d20_userreg1e-06_itemreg1e-05_lr01_CT20211203 gwas_int_data_for_pr_CT20211130 'test'
echo "Done running"
