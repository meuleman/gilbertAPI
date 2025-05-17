#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/job.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/job.out
pwd; hostname; date

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data/rank_based

#####Variables####
GENES_FILENAME=./gencode.v41.basic.annotation.parsed.promoterbp.extendeddomain.CT20221026.txt
GENES_CLEANED_FILENAME=./gencode.v41.basic.annotation.parsed.promoterbp.extendeddomain.CT20221027.cleaned.bed
GENES_CLEANED_C_FILENAME=./gencode.v41.basic.annotation.parsed.promoterbp.extendeddomain.CT20221027.cleaned.c1_11_12_14.bed
GENES_CLEANED_C_SORTED_FILENAME=./gencode.v41.basic.annotation.parsed.promoterbp.extendeddomain.CT20221027.cleaned.c1_11_12_14.sorted.bed

SORTED_SCALE_FILENAME=../Altius_25kb_hg38_maxMean_CT20220815.sorted.c1t3.bed
NEW_MAPPED_FILE=regions_25kb_CT20221027.proteinidmapped.bed

RANKED_REGIONS_OUT=./regions_25kb_CT20221109.STRINGsimilarityrankedregions.hdf5

NUM_CPUS=10
##################


##cleaning gene dataset
#python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_clean_gene_data_CT20221027.py $GENES_FILENAME $GENES_CLEANED_FILENAME
#cut -f 1,11,12,14 $GENES_CLEANED_FILENAME > $GENES_CLEANED_C_FILENAME
#sort-bed $GENES_CLEANED_C_FILENAME > $GENES_CLEANED_C_SORTED_FILENAME

##mapping features to sorted scale file
#bedmap --echo --echo-map-id $SORTED_SCALE_FILENAME $GENES_CLEANED_C_SORTED_FILENAME > $NEW_MAPPED_FILE

#echo "Generating each metric"

stringSim=()
#coexpression_human_transferred_max cooccurence database_human_transferred_max experiments_human_transferred_max fusion homology neighborhood_human_transferred_max textmining_human_transferred_max
for metric in coexpression_human_transferred_max cooccurence database_human_transferred_max experiments_human_transferred_max fusion homology neighborhood_human_transferred_max textmining_human_transferred_max
do
#   echo $metric
#   INFILE="../STRING/9606.protein.links.v11.5."$metric".all_associations.CT20220711.txt"
   OUTHDF5="./regions_25kb_CT20221031."$metric"STRINGsimilarity.hdf5"
   stringSim+=($OUTHDF5)
#   python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_STRING_similarity_CT20221031.py $INFILE $NEW_MAPPED_FILE $OUTHDF5 $NUM_CPUS
done


echo "Ranking regions"
python /home/ctrader/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_rank_regional_gene_interconnections_CT20221109.py "${stringSim[@]}" $RANKED_REGIONS_OUT

echo "Done"
