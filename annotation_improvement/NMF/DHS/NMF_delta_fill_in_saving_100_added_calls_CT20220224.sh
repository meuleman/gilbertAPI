#!/bin/bash
#SBATCH --job-name=fill_in
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50000
#SBATCH --error=/home/ctrader/Hold/NMF/DHS/slurm_output/job.%J.err
#SBATCH --output=/home/ctrader/Hold/NMF/DHS/slurm_output/job.%J.out
pwd; hostname; date

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/AG/AugmentedGenomics/recommendersystems/NMF/DHS/

echo "Running job array"
python NMF_delta_fill_in_saving_100_added_calls_CT20220224.py
echo "Done running"
