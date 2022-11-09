#!/bin/bash
#SBATCH --job-name=domain_extension
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=15000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/gencode_domain_extension.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/gencode_domain_extension.out
pwd; hostname; date

module load anaconda

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data

#####Variables####
IN_PATH=/home/ctrader/region_similarity/data/gencode.v41.basic.annotation.gtf
OUT_PATH=/home/ctrader/region_similarity/data/gencode.v41.basic.annotation.parsed.promoterbp.extendeddomain.CT20221026.txt

NUM_CPUS=15
##################

python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/gencode_domain_extension_CT20221026.py $IN_PATH $OUT_PATH $NUM_CPUS

echo "Done"
