#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=250000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/sfc_preprocessing_full.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/sfc_preprocessing_full.out
pwd; hostname; date
today=$(date "+%Y%m%d")

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm


indexed_dir=/home/ctrader/region_similarity/data/rank_based/indexed/hilbert_sfc/${today}/
output_dir=/home/ctrader/region_similarity/data/rank_based/hilbert/sfc/${today}/
mkdir ${output_dir}
mkdir ${indexed_dir}
cd ${output_dir}


for order in 4 5 6 7 8 9 10 11
do
  echo $order

  ####general variables####
  hilbert_dir=/home/ctrader/hilbert/regions/
  sorted_scale_file=${hilbert_dir}hilbert_order_${order}_windows.sorted.bed
  sorted_rownum4=hilbert_order_${order}_windows.sorted.rownum4.bed
  sorted_ranks5t38=hilbert_order_${order}_windows.sorted.rownum4.sfc.ranks5t38.bed

  #####cs variables####
  cs_datadir=/home/ctrader/hilbert/results/chromatin_states/20230622/
  cs_mapped_file=${cs_datadir}order_${order}_windows_annotation_sfc_scores.txt.gz
  cs_base_hdf5_file=hilbert_regions_order_${order}_CT${today}.cs.sfc.hdf5
  cs_percentile_file=hilbert_regions_order_${order}_CT${today}.cs.sfc.ranks.hdf5

  ####dhs variables####
  dhs_datadir=/home/ctrader/hilbert/results/DHS/20230622/
  dhs_mapped_file=${dhs_datadir}order_${order}_windows_annotation_sfc_scores.txt.gz
  dhs_base_hdf5_file=hilbert_regions_order_${order}_CT${today}.dhs.sfc.hdf5
  dhs_percentile_file=hilbert_regions_order_${order}_CT${today}.dhs.sfc.ranks.hdf5

  ####joined metrics variables####
  joined_metrics_file=hilbert_regions_order_${order}_CT${today}.dhs.cs.sfc.ranks.hdf5
  joined_metrics_txt=hilbert_regions_order_${order}_CT${today}.dhs.cs.sfc.ranks.txt
  top_region_thresh=21
  multi_factor_thresh=2
  search_space_file=hilbert_regions_order_${order}_CT${today}.dhs.cs.sfc.ranks.search_space.ft_${multi_factor_thresh}.rt_${top_region_thresh}.npz

  ####cs####
  # SFC values
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/sfc_to_hdf5_CT20230718.py $cs_mapped_file $cs_base_hdf5_file 'cs'
  # Ranking
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/rank_sfc_CT20230718.py $cs_base_hdf5_file $cs_percentile_file 'cs'

  ####dhs####
  # SFC values
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/sfc_to_hdf5_CT20230718.py $dhs_mapped_file $dhs_base_hdf5_file 'dhs'
  # Ranking
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/rank_sfc_CT20230718.py $dhs_base_hdf5_file $dhs_percentile_file 'dhs'

  ####Joining metrics####
  echo "Joining metrics"
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/joining_metrics_CT20230628.py $dhs_percentile_file $cs_percentile_file $sorted_scale_file $joined_metrics_file $joined_metrics_txt
  echo "Limiting search space"
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/search_space_CT20230718.py $dhs_base_hdf5_file $cs_base_hdf5_file $joined_metrics_file $top_region_thresh $multi_factor_thresh $search_space_file


  ####Adding Ranks to Tabix Files####
  awk '{print $0,"\t"NR-1}' $sorted_scale_file > $sorted_rownum4
  paste -d"\t" $sorted_rownum4 $joined_metrics_txt > $sorted_ranks5t38
  bgzip $sorted_ranks5t38 > $sorted_ranks5t38.bgz
  tabix -s 1 -b 2 -e 3 -f $sorted_ranks5t38.gz

  cp $sorted_ranks5t38.gz ${indexed_dir}
  cp $sorted_ranks5t38.gz.tbi ${indexed_dir}

done

echo "Done"
