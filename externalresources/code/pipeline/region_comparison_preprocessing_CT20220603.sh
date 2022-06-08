#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=30000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/job.%J.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/job.%J.out

#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/job.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/job.out
pwd; hostname; date

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data

#####Variables####
SCALE_FILENAME=Altius_100k_epilogos_hg19_WM20210518.bed
SORTED_SCALE_FILENAME=Altius_100k_epilogos_hg19_WM20210518.sorted.bed
#MAPPING_DIR="./tf_motif_bed/*"
MAPPING_DIR="./dhs_component_bed/*"
#NEW_MAPPED_FILE=regions_100k_CT20220602.TFmapped.bed
NEW_MAPPED_FILE=regions_100k_CT20220602.DHSmapped.bed
#OOE_NPY_FILE=regions_100k_CT20220602.TFoee.npy
OOE_NPY_FILE=regions_100k_CT20220602.DHSoee.npy
COSINE_NPZ_FILE=./regions_100k_CT20220602.DHScosine.npz
#COSINE_NPZ_FILE=./regions_100k_CT20220602.TFcosine.npz
TMP_FILE=./DHS_tmp_files
NUM_CPUS=8
##################

#sorting bed files
echo "Sorting"
sort-bed $SCALE_FILENAME > $SORTED_SCALE_FILENAME

#mapping features to sorted scale file
#initializing new file
cut -f 1-3 $SORTED_SCALE_FILENAME > $NEW_MAPPED_FILE

echo "Mapping"
for f in $MAPPING_DIR
do
    MAP_NAME=${f##*/}
    MAP_NAME=${MAP_NAME%.*}
    f_SORTED=$MAP_NAME.sorted.bed
    sort-bed $f > $f_SORTED
    sed -i "s/$/\t${MAP_NAME}/" $NEW_MAPPED_FILE
    TMP_FILE1=${SORTED_SCALE_FILENAME%%.*}.sorted.$MAP_NAME.mapped.bed
    TMP_FILE2=mapping_tmp.bed
    bedmap --echo --count $SORTED_SCALE_FILENAME $f_SORTED > $TMP_FILE1
    cut -d'|' -f 2 $TMP_FILE1 | paste $NEW_MAPPED_FILE - > $TMP_FILE2
    mv $TMP_FILE2 ./$NEW_MAPPED_FILE
    rm $TMP_FILE1 $f_SORTED
done

#Calculating OOE
echo "Calculating OOE"
python /home/ctrader/AG/AugmentedGenomics/externalresources/code/region_comparison_ooe_CT20220602.py $NEW_MAPPED_FILE ${MAPPING_DIR%/*} $OOE_NPY_FILE

#Cosine Similarity
echo "Calculating Cosine Similarity"
python /home/ctrader/AG/AugmentedGenomics/externalresources/code/region_comparison_cosine_sim_CT20220603.py $OOE_NPY_FILE $NUM_CPUS $COSINE_NPZ_FILE $TMP_FILE

echo "Done"
