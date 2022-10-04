#!/bin/bash
#SBATCH --job-name=fill_in
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120000
#SBATCH --error=/home/ctrader/Hold/NMF/DHS/slurm_output/job.%J.err
#SBATCH --output=/home/ctrader/Hold/NMF/DHS/slurm_output/job.%J.out
pwd; hostname; date

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/AG/AugmentedGenomics/recommendersystems/NMF/DHS/

echo "Running job array"
python NMF_delta_flip_by_unique_delta_vals_CT20220301.py 10 100
echo "Done running"
