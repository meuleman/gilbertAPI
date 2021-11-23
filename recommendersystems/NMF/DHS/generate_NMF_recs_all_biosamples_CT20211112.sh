#!/bin/bash
#SBATCH --job-name=recs_generation
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=50000
#SBATCH --error=/home/ctrader/Hold/DHS/slurm_output/job.%J.err
#SBATCH --output=/home/ctrader/Hold/DHS/slurm_output/job.%J.out
pwd; hostname; date

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/AG/AugmentedGenomics/recommendersystems/NMF/DHS/

echo "Running job array"
python generate_NMF_recs_all_biosamples_CT20211112.py
echo "Done running"
