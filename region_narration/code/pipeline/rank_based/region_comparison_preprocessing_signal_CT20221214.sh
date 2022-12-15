#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/signal_preprocessing.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/signal_preprocessing.out
pwd; hostname; date

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data/rank_based

#####Variables####
SORTED_SCALE_FILENAME=../Altius_25kb_hg38_maxMean_CT20220815.sorted.c1t3.bed

#DHS
#MAPPING_DIR="../dhs_component_bed/*"
#NEW_MAPPED_FILE=regions_25kb_CT20221207.DHSsignalmapped.bed
#Z_HDF5_FILE=regions_25kb_CT20221207.DHSsignalzscores.hdf5
#CONV_HDF5_FILE=./regions_25kb_CT20221214.DHSsignal.z.per.conv.sqrt.sigmoid.hdf5

#TF Motifs
#MAPPING_DIR="../tf_motif_bed_CT20221213/*"
#NEW_MAPPED_FILE=regions_25kb_CT20221213.TFsignalmapped.bed
Z_HDF5_FILE=regions_25kb_CT20221214.TFsignalzscores.hdf5
CONV_HDF5_FILE=./regions_25kb_CT20221214.TFsignal.z.per.conv.sqrt.sigmoid.hdf5
##################

#mapping features to sorted scale file
#initializing new file

#cut -f 1-3 $SORTED_SCALE_FILENAME > $NEW_MAPPED_FILE
#HEADER="chrm\tstart\tend"

#echo "Mapping"
#for f in $MAPPING_DIR
#do
#    MAP_NAME=${f##*/}
#    MAP_NAME=${MAP_NAME%.*}
#    HEADER+="\t$MAP_NAME"
#    f_SORTED=$MAP_NAME.sorted.bed
#    if [[ $MAPPING_DIR == *"dhs"* ]]; then
#        echo "DHS test worked"
#        sort-bed $f | awk -F '\t' -v OFS="\t" '{print $1, $2, $3, $9, $4}' | awk -F '\t' -v OFS="\t" '{gsub(" ", "", $4)}1' > $f_SORTED
#    else
#        echo "TF test worked"
#        sort-bed $f | awk -F '\t' -v OFS="\t" '{print $1, $2, $3, $4, $5}' | awk -F '\t' -v OFS="\t" '{gsub("/", "_", $4)}1' > $f_SORTED
#    fi
#    TMP_FILE2=$MAP_NAME.tmp.bed
#    bedmap --echo --count --mean --echo-map-score  $SORTED_SCALE_FILENAME $f_SORTED | cut -d'|' -f 2,3,4 | paste $NEW_MAPPED_FILE - > $TMP_FILE2
#    mv $TMP_FILE2 ./$NEW_MAPPED_FILE
#    rm $f_SORTED
#done

#sed -i "1i $HEADER" $NEW_MAPPED_FILE

##counts
#echo "Z Score"
#python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_signal_zscore_CT20221214.py  $NEW_MAPPED_FILE $Z_HDF5_FILE

#Sigmoid Euclidean Similarity
echo "Calculating Similarity"
python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_zscore_sim_CT20221214.py $Z_HDF5_FILE $CONV_HDF5_FILE

echo "Done"
