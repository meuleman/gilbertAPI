#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/preprocessing_full.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/preprocessing_full.out
pwd; hostname; date
today=$(date "+%Y%m%d")

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data/rank_based


for scale in 1000 2000 5000 10000 15000 25000 50000 75000 100000
do
  echo $scale

  ####general variables####
  chrm_length_file=hg38_chromosome_lengths_CT20230327.tsv
  scale_file=hg38.${scale/$1000}kb.regions.CT${today}.bed
  sorted_scale_file=hg38.${scale/$1000}kb.regions.CT${today}.sorted.bed
  sorted_rownum4=hg38.${scale/$1000}kb.regions.CT${today}.sorted.rownum4.bed
  bw=1000

  #####epilogos variables####

  epi_mapping_dir="epilogos_component_bed_CT20230203/*"
  epi_new_mapped_file=regions_${scale/$1000}kb_CT${today}.epilogos.summapped.bed
  epi_base_hdf5_file=regions_${scale/$1000}kb_CT${today}.epilogos.statesums.hdf5
  epi_percentile_file=regions_${scale/$1000}kb_CT${today}.epilogos.statesum_rank_percentile.hdf5

  ####dhs variables####

  dhs_mapping_dir="DHS_NMF_factor_norm_bed_CT20230210/*"
  dhs_new_mapped_file=regions_${scale/$1000}kb_CT${today}.DHS.factornormNMF.summapped.bed
  dhs_base_hdf5_file=regions_${scale/$1000}kb_CT${today}.DHS.factornormNMF.factorsums.hdf5
  dhs_percentile_file=regions_${scale/$1000}kb_CT${today}.DHS.factornormNMF.sum_enrich_corank_percentile.hdf5

  ####joined metrics variables####

  joined_metrics_file=regions_${scale/$1000}kb_CT${today}.DHS.factornormNMF.sum_enrich_corank.epilogos.statesum_rank.percentiles.hdf5
  search_space_thresh=0.95
  search_space_file=regions_${scale/$1000}kb_CT${today}.DHS.factornormNMF.sum_enrich_corank.epilogos.statesum_rank.percentiles.limited_search_space.npz




  ####Create regions file####
  #create regions of size:scale and windowed by size:bw across genome
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/creating_regions_by_scale_CT20230327.py $scale $bw $chrm_length_file $scale_file
  sort-bed $scale_file > $sorted_scale_file

  awk '{print $0,"\t"NR-1}' $sorted_scale_file > $sorted_rownum4
  bgzip $sorted_rownum4 > $sorted_rownum4.bgz
  tabix -s 1 -b 2 -e 3 $sorted_rownum4.gz




  ####epilogos mapping####

  #mapping features to sorted scale file
  #initializing new file

  cut -f 1-3 $sorted_scale_file > $epi_new_mapped_file
  header="chrm\tstart\tend"

  echo "Epilogos Mapping"
  for f in $epi_mapping_dir
  do
      map_name=${f##*/}
      map_name=${map_name%.*}
      header+="\t$map_name"
      f_sorted=$map_name.sorted.bed
      sort-bed $f > $f_sorted
      tmp_file=$map_name.tmp.bed
      bedmap --echo --sum $sorted_scale_file $f_sorted | cut -d'|' -f 2 | paste $epi_new_mapped_file - > $tmp_file
      mv $tmp_file ./$epi_new_mapped_file
      rm $f_sorted
  done

  sed -i "1i $header" $epi_new_mapped_file

  # Means
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_epilogos_statesum_CT20230209.py $epi_new_mapped_file $epi_base_hdf5_file

  # Ranking
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_epilogos_sum_rank_CT20230214.py $epi_base_hdf5_file $epi_percentile_file




  ####dhs mapping####

  #mapping features to sorted scale file
  #initializing new file

  cut -f 1-3 $sorted_scale_file > $dhs_new_mapped_file
  header="chrm\tstart\tend"

  echo "DHS Mapping"
  for f in $dhs_mapping_dir
  do
      map_name=${f##*/}
      map_name=${map_name%.*}
      header+="\t$map_name"
      f_sorted=$map_name.sorted.bed
      sort-bed $f > $f_sorted
      tmp_file=$map_name.tmp.bed
      bedmap --echo --sum $sorted_scale_file $f_sorted | cut -d'|' -f 2 | paste $dhs_new_mapped_file - > $tmp_file
      mv $tmp_file ./$dhs_new_mapped_file
      rm $f_sorted
  done

  sed -i "1i $header" $dhs_new_mapped_file

  # Sum processing
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_factorsum_CT20230210.py $dhs_new_mapped_file $dhs_base_hdf5_file

  # Ranking
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_DHS_NMF_sum_enrich_corank_CT20230214.py $dhs_base_hdf5_file $dhs_percentile_file





  ####Joining metrics####
  echo "Joining metrics"
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/joining_metrics_CT20230327.py $dhs_percentile_file $epi_percentile_file $sorted_scale_file $joined_metrics_file

  echo "Limiting search space"
  python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/limiting_search_space_CT20230327.py $joined_metrics_file $search_space_thresh $search_space_file

done

echo "Done"
