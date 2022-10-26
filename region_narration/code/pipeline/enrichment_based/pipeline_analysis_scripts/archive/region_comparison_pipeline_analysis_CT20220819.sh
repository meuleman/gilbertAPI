#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/pipeline_analysis.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/pipeline_analysis.out
pwd; hostname; date

module load anaconda
#module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data

#####Variables####
SORTED_SCALE_FILENAME=Altius_25kb_hg38_maxMean_CT20220815.sorted.c1t3.bed
SCALE_GZ_FILENAME=Altius_25kb_hg38_maxMean_CT20220815.sorted.c1t3.rownum4.bed.gz

NUM_CPUS=1

DHS_OOE_HDF5_FILE=regions_25kb_CT20220815.DHSoee.hdf5
TF_OOE_HDF5_FILE=regions_25kb_CT20220815.TFoee.hdf5
EPILOGOS_OOE_HDF5_FILE=regions_25kb_CT20220815.epilogosoee.hdf5

DHS_MAPPED_FILE=regions_25kb_CT20220815.DHSmapped.bed
TF_MAPPED_FILE=regions_25kb_CT20220815.TFmapped.bed
EPILOGOS_MAPPED_FILE=regions_25kb_CT20220815.epilogosmapped.bed

DHS_COSINE_HDF5_FILE=./regions_25kb_CT20220815.DHScosine.hdf5
TF_COSINE_HDF5_FILE=./regions_25kb_CT20220815.TFcosine.hdf5
EPILOGOS_COSINE_HDF5_FILE=./regions_25kb_CT20220815.epilogoscosine.hdf5
COEX_HDF5_FILE=./regions_25kb_CT20220815.coexpression_human_transferred_maxSTRINGsimilarity.hdf5
COOCCUR_HDF5_FILE=./regions_25kb_CT20220815.cooccurenceSTRINGsimilarity.hdf5
DATABASE_HDF5_FILE=./regions_25kb_CT20220815.database_human_transferred_maxSTRINGsimilarity.hdf5
EXPER_HDF5_FILE=./regions_25kb_CT20220815.experiments_human_transferred_maxSTRINGsimilarity.hdf5
FUSION_HDF5_FILE=./regions_25kb_CT20220815.fusionSTRINGsimilarity.hdf5
HOMOLOGY_HDF5_FILE=./regions_25kb_CT20220815.homologySTRINGsimilarity.hdf5
NEIGH_HDF5_FILE=./regions_25kb_CT20220815.neighborhood_human_transferred_maxSTRINGsimilarity.hdf5
TEXTMIN_HDF5_FILE=./regions_25kb_CT20220815.textmining_human_transferred_maxSTRINGsimilarity.hdf5

DHSFACT2DESCR_FILE=./DHS_descriptor_dict_CT20220729.pickle
TFFACT2DESCR_FILE=./tf_motif_descriptor_dict_CT20220627.pickle
EPIFACTOR2DESCR_FILE=./epilogos_descriptor_dict_CT20220627.pickle

GENEID2NAME_PATH=./gene_proteinID2genename_dict_CT20220712.pickle

GENE_SET_ORDER_PATH=./regions_25kb_CT20220817.genesethallmarkc1c2order.cosine.average.npz
GENE_SET_ALL=./regions_25kb_CT20220817.genesethallmarkc1c2freqooe.hdf5
GENE_SET_SIG=./regions_25kb_CT20220817.genesethallmarkc1c2freqooe.sig.hdf5


OUT_BED_PATH=/home/ctrader/public_html/
##################
GENOMIC_LOCATION=$1

SENTENCE_OUT=./sentence_summaries/regions_25kb_${GENOMIC_LOCATION/:/_}_sentence_summaries_CT20220818.txt
UPSET_LIMITING_PER=95
##################

echo $GENOMIC_LOCATION

#find similar regions
echo "Generating Sentences"
time python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/pipeline_analysis_scripts/region_comparison_pipeline_analysis_sentences_CT20220818.py $SCALE_GZ_FILENAME \
$GENOMIC_LOCATION $DHS_COSINE_HDF5_FILE $TF_COSINE_HDF5_FILE $EPILOGOS_COSINE_HDF5_FILE $COEX_HDF5_FILE $COOCCUR_HDF5_FILE $DATABASE_HDF5_FILE $EXPER_HDF5_FILE $FUSION_HDF5_FILE \
$HOMOLOGY_HDF5_FILE $NEIGH_HDF5_FILE $TEXTMIN_HDF5_FILE $SORTED_SCALE_FILENAME $DHS_OOE_HDF5_FILE $TF_OOE_HDF5_FILE $EPILOGOS_OOE_HDF5_FILE $DHS_MAPPED_FILE $TF_MAPPED_FILE \
$EPILOGOS_MAPPED_FILE $DHSFACT2DESCR_FILE $TFFACT2DESCR_FILE $EPIFACTOR2DESCR_FILE $GENEID2NAME_PATH $OUT_BED_PATH $SENTENCE_OUT

echo "Making bar plots"
time python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/pipeline_analysis_scripts/region_comparison_pipeline_analysis_bar_CT20220818.py $SCALE_GZ_FILENAME \
$GENOMIC_LOCATION $DHS_COSINE_HDF5_FILE $TF_COSINE_HDF5_FILE $EPILOGOS_COSINE_HDF5_FILE $COEX_HDF5_FILE $COOCCUR_HDF5_FILE $DATABASE_HDF5_FILE $EXPER_HDF5_FILE $FUSION_HDF5_FILE \
$HOMOLOGY_HDF5_FILE $NEIGH_HDF5_FILE $TEXTMIN_HDF5_FILE $SORTED_SCALE_FILENAME $DHS_OOE_HDF5_FILE $TF_OOE_HDF5_FILE $EPILOGOS_OOE_HDF5_FILE $DHS_MAPPED_FILE $TF_MAPPED_FILE \
$EPILOGOS_MAPPED_FILE $DHSFACT2DESCR_FILE $TFFACT2DESCR_FILE $EPIFACTOR2DESCR_FILE $GENEID2NAME_PATH $OUT_BED_PATH

echo "Making upset plots"
time python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/pipeline_analysis_scripts/region_comparison_pipeline_analysis_upset_CT20220819.py $SCALE_GZ_FILENAME \
$GENOMIC_LOCATION $DHS_COSINE_HDF5_FILE $TF_COSINE_HDF5_FILE $EPILOGOS_COSINE_HDF5_FILE $COEX_HDF5_FILE $COOCCUR_HDF5_FILE $DATABASE_HDF5_FILE $EXPER_HDF5_FILE $FUSION_HDF5_FILE \
$HOMOLOGY_HDF5_FILE $NEIGH_HDF5_FILE $TEXTMIN_HDF5_FILE $SORTED_SCALE_FILENAME $DHS_OOE_HDF5_FILE $TF_OOE_HDF5_FILE $EPILOGOS_OOE_HDF5_FILE $DHS_MAPPED_FILE $TF_MAPPED_FILE \
$EPILOGOS_MAPPED_FILE $DHSFACT2DESCR_FILE $TFFACT2DESCR_FILE $EPIFACTOR2DESCR_FILE $GENEID2NAME_PATH $OUT_BED_PATH $UPSET_LIMITING_PER

echo "Making spectrum plots"
time python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/pipeline_analysis_scripts/region_comparison_pipeline_analysis_spectrum_CT20220819.py $SCALE_GZ_FILENAME \
$GENOMIC_LOCATION $SORTED_SCALE_FILENAME $GENE_SET_ORDER_PATH $GENE_SET_SIG $GENE_SET_ALL

echo "Done"
