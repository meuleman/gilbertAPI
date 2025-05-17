#!/bin/bash
#SBATCH --job-name=corrs_CLT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40000
#SBATCH --error=/home/ctrader/megamap/slurm_out/clt.err
#SBATCH --output=/home/ctrader/megamap/slurm_out/clt.out
pwd; hostname; date

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

#####Variables####
HALLMARK=/home/ctrader/region_similarity/data/gene2hallmarkgeneset_dict_CT20220718.pickle
C1=/home/ctrader/region_similarity/data/gene2c1geneset_dict_CT20220725.pickle
C2=/home/ctrader/region_similarity/data/gene2c2geneset_dict_CT20220725.pickle
#CORRS=/net/seq/data2/projects/megamap/rnaseq_analysis_v4/metric_thresh_cossim/deseq_stat/effect_thresh_2.0/min_num_effect_3/corrs.feather
CORRS=/net/seq/data2/projects/megamap/rnaseq_analysis_v4/metric_pvals_cossim/deseq_lfc_shrunk_prior_0.25/effect_thresh_2.0/min_num_effect_2/corrs.feather

#OUT=/home/ctrader/megamap/data/geneset_significance.CLT.v4.cossim.effect_thresh_2.min_num_effect_3.hdf5
OUT=/home/ctrader/megamap/data/geneset_significance.CLT.v4.pvalcossim.shrunkprior0_25.effect_thresh_2.min_num_effect_2.hdf5
##################
python ./gene_set_significance_clt_CT20221012.py $CORRS $HALLMARK $C1 $C2 $OUT
