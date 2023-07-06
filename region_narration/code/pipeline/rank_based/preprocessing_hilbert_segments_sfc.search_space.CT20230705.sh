#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=250000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/sfc_preprocessing_subset.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/sfc_preprocessing_subset.out
pwd; hostname; date
today=$(date "+%Y%m%d")

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data/rank_based/hilbert/sfc/


for order in 4 5 6 7 8 9
do
  echo $order

  ####joined metrics variables####

  joined_metrics_file=hilbert_regions_order_${order}_CT20230630.dhs.cs.sfc.ranks.hdf5
  search_space_thresh=1000
  multi_factor_thresh=2
  search_space_file=hilbert_regions_order_${order}_CT${today}.dhs.cs.sfc.ranks.limited_search_space.top_${search_space_thresh}_per_factor.gte_${multi_factor_thresh}_ranked_factors.npz


  echo "Limiting search space"
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/limiting_search_space_with_expansion_CT20230705.py $joined_metrics_file $search_space_thresh $multi_factor_thresh $search_space_file


done

echo "Done"
