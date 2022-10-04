#!/bin/bash
#SBATCH --job-name=fill_in
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=80000
#SBATCH --error=/home/ctrader/Hold/DHS/slurm_output/job.%J.err
#SBATCH --output=/home/ctrader/Hold/DHS/slurm_output/job.%J.out
pwd; hostname; date

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/AG/AugmentedGenomics/recommendersystems/NMF/DHS/

echo "Running job array"
python NMF_delta_fill_in_one_at_a_time_local_NMF_CT20220209.py
echo "Done running"
