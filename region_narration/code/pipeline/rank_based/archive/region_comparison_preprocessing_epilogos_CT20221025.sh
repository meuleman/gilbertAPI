#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/preprocessing_epilogos.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/preprocessing_epilogos.out
pwd; hostname; date

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data/count_based

#####Variables####
SORTED_SCALE_FILENAME=../Altius_25kb_hg38_maxMean_CT20220815.sorted.c1t3.bed

MAPPING_DIR="../epilogos_component_bed_CT20220815/*"
NEW_MAPPED_FILE=regions_25kb_CT20221025.epilogosmapped.bed
BASE_HDF5_FILE=regions_25kb_CT20221025.epilogosmeanmax.hdf5
SIGMOIDEUC_HDF5_FILE=./regions_25kb_CT20221025.epilogossigmoideuc.hdf5

NUM_CPUS=1
CS=True
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
    bedmap --echo --mean --max $SORTED_SCALE_FILENAME $f_SORTED | cut -d'|' -f 2,3 | paste $NEW_MAPPED_FILE - > $TMP_FILE2
    mv $TMP_FILE2 ./$NEW_MAPPED_FILE
    rm $f_SORTED
done

sed -i "1i $HEADER" $NEW_MAPPED_FILE

#Calculating OOE
echo "Calculating OOE"
python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/count_based/region_comparison_meanmax_epilogos_CT20221025.py $NEW_MAPPED_FILE $BASE_HDF5_FILE

#Sigmoid Euclidean Similarity
echo "Calculating Sigmoid Euclidean Similarity"
python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/count_based/region_comparison_sigmoid_eucdist_CT20221025.py $BASE_HDF5_FILE $NUM_CPUS $SIGMOIDEUC_HDF5_FILE $CS

echo "Done"
