#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/job.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/job.out
pwd; hostname; date

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data

#####Variables####
GENES_FILENAME=./mart_export_protein_coding_genes_unique_with_proteinID_CT20220701.txt
GENES_CLEANED_FILENAME=./mart_export_protein_coding_genes_unique_with_proteinID_CT20220712.cleaned.bed
GENES_CLEANED_SORTED_FILENAME=./mart_export_protein_coding_genes_unique_with_proteinID_CT20220712.cleaned.sorted.bed
GENES_CLEANED_SORTED_C1236_FILENAME=./mart_export_protein_coding_genes_unique_with_proteinID_CT20220712.cleaned.sorted.c1236.bed

SORTED_SCALE_FILENAME=Altius_25kb_hg38_maxMean_CT20220613.sorted.c1t3.bed
NEW_MAPPED_FILE=regions_25kb_CT20220712.geneidmapped.bed

NUM_CPUS=10
##################


#cleaning gene dataset
python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_clean_gene_data_CT20220705.py $GENES_FILENAME $GENES_CLEANED_FILENAME
sed -i 's/^/chr/' $GENES_CLEANED_FILENAME
sort-bed $GENES_CLEANED_FILENAME > $GENES_CLEANED_SORTED_FILENAME

cut -f 1-3,6 $GENES_CLEANED_SORTED_FILENAME > $GENES_CLEANED_SORTED_C1236_FILENAME

#mapping features to sorted scale file
bedmap --echo --echo-map-id $SORTED_SCALE_FILENAME $GENES_CLEANED_SORTED_C1236_FILENAME > $NEW_MAPPED_FILE

echo "Generating each metric"

for metric in coexpression_human_transferred_max cooccurence database_human_transferred_max experiments_human_transferred_max fusion homology neighborhood_human_transferred_max textmining_human_transferred_max
do
   echo $metric
   INFILE="./STRING/9606.protein.links.v11.5."$metric".all_associations.CT20220711.txt"
   OUTHDF5="./regions_25kb_CT20220712."$metric"STRINGsimilarity.hdf5"
   #OUTNPZ="./regions_25kb_CT20220712."$metric"STRINGoriginalgenes.npz"
   python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_STRING_similarity_CT20220712.py $INFILE $NEW_MAPPED_FILE $OUTHDF5 $NUM_CPUS
done

echo "Done"
