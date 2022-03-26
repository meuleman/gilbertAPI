#!/bin/bash
#SBATCH --job-name=delta_flip
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100000
#SBATCH --error=/home/ctrader/Hold/NMF/DHS/slurm_output/job.%J.err
#SBATCH --output=/home/ctrader/Hold/NMF/DHS/slurm_output/job.%J.out
pwd; hostname; date

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/AG/AugmentedGenomics/recommendersystems/NMF/DHS/

echo "Running job array"
python NMF_delta_threshold_flip_CT20220314.py 1 -0.99 0.99
echo "Done running"
