#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/preprocessing_epilogos.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/preprocessing_epilogos.out
pwd; hostname; date

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data/rank_based

#####Variables####
SORTED_SCALE_FILENAME=hg38.25kb.maxMean.CT20230118.sorted.bed

MAPPING_DIR="DHS_NMF_factor_norm_bed_CT20230210/*"
NEW_MAPPED_FILE=regions_25kb_CT20230210.DHS.factornormNMF.summapped.bed
BASE_HDF5_FILE=regions_25kb_CT20230210.DHS.factornormNMF.factorsums.hdf5
PERCENTILE_FILE=regions_25kb_CT20230214.DHS.factornormNMF.sum_enrich_corank_percentile.hdf5
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
    bedmap --echo --sum $SORTED_SCALE_FILENAME $f_SORTED | cut -d'|' -f 2 | paste $NEW_MAPPED_FILE - > $TMP_FILE2
    mv $TMP_FILE2 ./$NEW_MAPPED_FILE
    rm $f_SORTED
done

sed -i "1i $HEADER" $NEW_MAPPED_FILE

# Sum processing
python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_factorsum_CT20230210.py $NEW_MAPPED_FILE $BASE_HDF5_FILE

# Ranking
python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_DHS_NMF_sum_enrich_corank_CT20230214.py $BASE_HDF5_FILE $PERCENTILE_FILE

echo "Done"
