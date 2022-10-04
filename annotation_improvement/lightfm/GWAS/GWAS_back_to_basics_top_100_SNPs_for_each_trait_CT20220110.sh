#!/bin/bash
#SBATCH --job-name=GWAS_top_100_SNPs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=15000
#SBATCH --error=/home/ctrader/Hold/GWAS/pickled_subthreshold_data/slurm_output/job.%J.err
#SBATCH --output=/home/ctrader/Hold/GWAS/pickled_subthreshold_data/slurm_output/job.%J.out
#SBATCH --array=0-3810:5
pwd; hostname; date

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/AG/AugmentedGenomics/recommendersystems/lightfm/GWAS/

echo "Running job array"
python GWAS_back_to_basics_top_100_SNPs_for_each_trait_CT20220110.py $SLURM_ARRAY_TASK_ID
echo "Done running"
