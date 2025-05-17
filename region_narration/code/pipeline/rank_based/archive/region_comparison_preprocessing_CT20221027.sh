#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/preprocessing.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/preprocessing.out
pwd; hostname; date

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data/count_based

#####Variables####
SORTED_SCALE_FILENAME=../Altius_25kb_hg38_maxMean_CT20220815.sorted.c1t3.bed

MAPPING_DIR="../tf_motif_bed/*"
NEW_MAPPED_FILE=regions_25kb_CT20221027.TFmapped.bed
COUNTS_HDF5_FILE=regions_25kb_CT20221027.TFcounts.hdf5
SIGMOIDEUC_HDF5_FILE=./regions_25kb_CT20221027.TFsigmoideuc.hdf5

#MAPPING_DIR="../dhs_component_bed/*"
#NEW_MAPPED_FILE=regions_25kb_CT20221027.DHSmapped.bed
#COUNTS_HDF5_FILE=regions_25kb_CT20221027.DHScounts.hdf5
#SIGMOIDEUC_HDF5_FILE=./regions_25kb_CT20221027.DHSsigmoideuc.hdf5

NUM_CPUS=1
##################

#mapping features to sorted scale file
#initializing new file
cut -f 1-3 $SORTED_SCALE_FILENAME > $NEW_MAPPED_FILE
HEADER="chrm\tstart\tend"

echo "Mapping"
for f in $MAPPING_DIR
do
    MAP_NAME=${f##*/}
    MAP_NAME=${MAP_NAME%.*}
    HEADER+="\t$MAP_NAME"
    f_SORTED=$MAP_NAME.sorted.bed
    sort-bed $f > $f_SORTED
    TMP_FILE2=$MAP_NAME.tmp.bed
    bedmap --echo --count $SORTED_SCALE_FILENAME $f_SORTED | cut -d'|' -f 2 | paste $NEW_MAPPED_FILE - > $TMP_FILE2
    mv $TMP_FILE2 ./$NEW_MAPPED_FILE
    rm $f_SORTED
done

sed -i "1i $HEADER" $NEW_MAPPED_FILE

#counts
echo "Counts"
python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/count_based/region_comparison_counts_CT20221025.py $NEW_MAPPED_FILE $COUNTS_HDF5_FILE

#Sigmoid Euclidean Similarity
echo "Calculating Similarity"
python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/count_based/region_comparison_sigmoid_eucdist_CT20221027.py $COUNTS_HDF5_FILE $NUM_CPUS $SIGMOIDEUC_HDF5_FILE

echo "Done"
