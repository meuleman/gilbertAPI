#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/gene_set_ooe.b.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/gene_set_ooe.b.out
pwd; hostname; date

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data

#####Variables####
EPILOGOS_COSINE_HDF5_FILE=./regions_25kb_CT20220815.epilogoscosine.hdf5
TF_COSINE_HDF5_FILE=./regions_25kb_CT20220815.TFcosine.hdf5
DHS_COSINE_HDF5_FILE=./regions_25kb_CT20220815.DHScosine.hdf5
COEX_HDF5_FILE=./regions_25kb_CT20220815.coexpression_human_transferred_maxSTRINGsimilarity.hdf5
COOCCUR_HDF5_FILE=./regions_25kb_CT20220815.cooccurenceSTRINGsimilarity.hdf5
DATABASE_HDF5_FILE=./regions_25kb_CT20220815.database_human_transferred_maxSTRINGsimilarity.hdf5
EXPER_HDF5_FILE=./regions_25kb_CT20220815.experiments_human_transferred_maxSTRINGsimilarity.hdf5
FUSION_HDF5_FILE=./regions_25kb_CT20220815.fusionSTRINGsimilarity.hdf5
HOMOLOGY_HDF5_FILE=./regions_25kb_CT20220815.homologySTRINGsimilarity.hdf5
NEIGH_HDF5_FILE=./regions_25kb_CT20220815.neighborhood_human_transferred_maxSTRINGsimilarity.hdf5
TEXTMIN_HDF5_FILE=./regions_25kb_CT20220815.textmining_human_transferred_maxSTRINGsimilarity.hdf5

SORTED_SCALE_FILENAME=Altius_25kb_hg38_maxMean_CT20220815.sorted.c1t3.bed
GENES_PRECLEANED=./mart_export_protein_coding_genes_unique_CT20220628.txt
GENES_BED=./mart_export_protein_coding_genes_unique_cleaned_CT20220817.bed
GENES_BED_SORTED=./mart_export_protein_coding_genes_unique_cleaned_CT20220817.sorted.bed
GENES_MAPPED_BED=./regions_25kb_CT20220817.cleanedgenesmapped.bed

GENE2GENESET_FILE=./gene2hallmarkc1c2genesets_dict_CT20220729.pickle

GENE_SET_OOE=./regions_25kb_CT20220817.genesethallmarkc1c2freqooe.hdf5

GENE_SET_OOE_SIG=./regions_25kb_CT20220817.genesethallmarkc1c2freqooe.sig.hdf5
#GENESET_ORDER_THROUGH_CLUSTERING_NPY=./regions_25kb_CT20220817.genesethallmarkc1c2order.cosine.average.npz
GENESET_ORDER_THROUGH_CLUSTERING_NPY=./regions_25kb_CT20220907.genesethallmarkc1c2order.significanceindependent.cosine.average.npz

NUM_CPUS=8

##################
#cleaning the gene file
#echo "Cleaning"
#python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_clean_gene_data_CT20220725.py $GENES_PRECLEANED $GENES_BED
##sed -i 's/^/chr/' $GENES_BED

#echo "Sorting"
#sort-bed $GENES_BED > $GENES_BED_SORTED

#echo "Mapping"
#bedmap --echo --echo-map $SORTED_SCALE_FILENAME $GENES_BED_SORTED > $GENES_MAPPED_BED

#Calculating OOE
#echo "Calculating OOE for Gene Sets"
#python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_gene_set_ooe_CT20220816.py $DHS_COSINE_HDF5_FILE $TF_COSINE_HDF5_FILE \
#$EPILOGOS_COSINE_HDF5_FILE $COEX_HDF5_FILE $COOCCUR_HDF5_FILE $DATABASE_HDF5_FILE $EXPER_HDF5_FILE $FUSION_HDF5_FILE $HOMOLOGY_HDF5_FILE $NEIGH_HDF5_FILE $TEXTMIN_HDF5_FILE \
#$GENES_MAPPED_BED $GENE2GENESET_FILE $GENE_SET_OOE $NUM_CPUS

#echo "Finding Significant Terms"
#python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_find_significant_gene_sets_CT20220817.py $GENE_SET_OOE $GENE_SET_OOE_SIG $NUM_CPUS

#echo "Finding Order of Terms through Clustering"
#python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_order_gene_sets_CT20220721.py $GENE_SET_OOE_SIG $GENESET_ORDER_THROUGH_CLUSTERING_NPY

echo "Finding Order of Terms through Clustering"
python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_order_gene_sets_significance_independent_CT20220907.py $GENE_SET_OOE $GENESET_ORDER_THROUGH_CLUSTERING_NPY

echo "Done"
