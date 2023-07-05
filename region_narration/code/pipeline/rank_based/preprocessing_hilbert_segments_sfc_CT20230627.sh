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

cd /home/ctrader/region_similarity/data/rank_based/hilbert/sfc/


for order in 4 5 6 7 8 9
do
  echo $order

  ####general variables####
  sorted_scale_file=hilbert_order_${order}_windows.sorted.bed
  sorted_rownum4=hilbert_order_${order}_windows.sorted.rownum4.bed
  sorted_ranks5t38=hilbert_order_${order}_windows.sorted.rownum4.sfc.ranks5t38.bed

  #####cs variables####
  cs_datadir=/home/ctrader/hilbert/results/chromatin_states/JQ_20230601/
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
  search_space_thresh=1000
  search_space_file=hilbert_regions_order_${order}_CT${today}.dhs.cs.sfc.ranks.limited_search_space_top_${search_space_thresh}_no_overlap.npz




  ####cs####

  # Means
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/sfc_to_hdf5_CT20230628.py $cs_mapped_file $cs_base_hdf5_file 'cs'

  # Ranking
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/rank_sfc_CT20230628.py $cs_base_hdf5_file $cs_percentile_file 'cs'




  ####dhs####

  # Sum processing
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/sfc_to_hdf5_CT20230628.py $dhs_mapped_file $dhs_base_hdf5_file 'dhs'

  # Ranking
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/rank_sfc_CT20230628.py $dhs_base_hdf5_file $dhs_percentile_file 'dhs'





  ####Joining metrics####
  echo "Joining metrics"
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/joining_metrics_CT20230628.py $dhs_percentile_file $cs_percentile_file $sorted_scale_file $joined_metrics_file $joined_metrics_txt

  echo "Limiting search space"
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/limiting_search_space_CT20230628.py $joined_metrics_file $search_space_thresh $search_space_file


  ####Adding Ranks to Tabix Files####
  awk '{print $0,"\t"NR-1}' $sorted_scale_file > $sorted_rownum4
  paste -d"\t" $sorted_rownum4 $joined_metrics_txt > $sorted_ranks5t38
  bgzip $sorted_ranks5t38 > $sorted_ranks5t38.bgz
  tabix -s 1 -b 2 -e 3 -f $sorted_ranks5t38.gz

  cp $sorted_ranks5t38.gz ../../indexed/hilbert_sfc/
  cp $sorted_ranks5t38.gz.tbi ../../indexed/hilbert_sfc/

done

echo "Done"
