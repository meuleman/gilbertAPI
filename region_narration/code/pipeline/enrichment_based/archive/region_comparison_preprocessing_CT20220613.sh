#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=20000
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
SORTED_SCALE_FILENAME=Altius_25kb_hg38_maxMean_CT20220613.sorted.c1t3.bed

#MAPPING_DIR="./tf_motif_bed/*"
#NEW_MAPPED_FILE=regions_25kb_CT20220613.TFmapped.bed
#OOE_NPY_FILE=regions_25kb_CT20220613.TFoee.npy
#COSINE_HDF5_FILE=./regions_25kb_CT20220613.TFcosine.hdf5

MAPPING_DIR="./dhs_component_bed/*"
NEW_MAPPED_FILE=regions_25kb_CT20220613.DHSmapped.bed
OOE_NPY_FILE=regions_25kb_CT20220613.DHSoee.npy
COSINE_HDF5_FILE=./regions_25kb_CT20220613.DHScosine.hdf5

NUM_CPUS=8
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
    #sed -i "s/$/\t${MAP_NAME}/" $NEW_MAPPED_FILE
    #TMP_FILE1=${SORTED_SCALE_FILENAME%%.*}.sorted.$MAP_NAME.mapped.bed
    TMP_FILE2=$MAP_NAME.tmp.bed
    bedmap --echo --count $SORTED_SCALE_FILENAME $f_SORTED | cut -d'|' -f 2 | paste $NEW_MAPPED_FILE - > $TMP_FILE2
    #cut -d'|' -f 2 $TMP_FILE1 | paste $NEW_MAPPED_FILE - > $TMP_FILE2
    mv $TMP_FILE2 ./$NEW_MAPPED_FILE
    rm $f_SORTED
done

sed -i "1i $HEADER" $NEW_MAPPED_FILE

#Calculating OOE
echo "Calculating OOE"
python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_ooe_CT20220608.py $NEW_MAPPED_FILE ${MAPPING_DIR%/*} $OOE_NPY_FILE

#Cosine Similarity
echo "Calculating Cosine Similarity"
python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_cosine_sim_CT20220610.py $OOE_NPY_FILE $NUM_CPUS $COSINE_HDF5_FILE

echo "Done"
