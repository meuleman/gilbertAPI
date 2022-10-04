#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/score_comb_job.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/score_comb_job.out
pwd; hostname; date

module load anaconda
#module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data

#####Variables####
SORTED_SCALE_FILENAME=Altius_25kb_hg38_maxMean_CT20220613.sorted.c1t3.bed
SCALE_GZ_FILENAME=Altius_25kb_hg38_maxMean_CT20220613.sorted.c1t3.rownum4.bed.gz

NUM_CPUS=1

DHS_OOE_HDF5_FILE=regions_25kb_CT20220623.DHSoee.hdf5
TF_OOE_HDF5_FILE=regions_25kb_CT20220623.TFoee.hdf5
EPILOGOS_OOE_HDF5_FILE=regions_25kb_CT20220623.epilogosoee.hdf5

DHS_MAPPED_FILE=regions_25kb_CT20220623.DHSmapped.bed
TF_MAPPED_FILE=regions_25kb_CT20220623.TFmapped.bed
EPILOGOS_MAPPED_FILE=regions_25kb_CT20220623.epilogosmapped.bed

DHS_COSINE_HDF5_FILE=./regions_25kb_CT20220623.DHScosine.hdf5
TF_COSINE_HDF5_FILE=./regions_25kb_CT20220623.TFcosine.hdf5
EPILOGOS_COSINE_HDF5_FILE=./regions_25kb_CT20220623.epilogoscosine.hdf5
COEX_HDF5_FILE=./regions_25kb_CT20220706.coexpression_human_transferred_maxSTRINGsimilarity.hdf5
COOCCUR_HDF5_FILE=./regions_25kb_CT20220706.cooccurenceSTRINGsimilarity.hdf5
DATABASE_HDF5_FILE=./regions_25kb_CT20220706.databaseSTRINGsimilarity.hdf5
EXPER_HDF5_FILE=./regions_25kb_CT20220706.experiments_human_transferred_maxSTRINGsimilarity.hdf5
FUSION_HDF5_FILE=./regions_25kb_CT20220706.fusionSTRINGsimilarity.hdf5
HOMOLOGY_HDF5_FILE=./regions_25kb_CT20220706.homologySTRINGsimilarity.hdf5
NEIGH_HDF5_FILE=./regions_25kb_CT20220706.neighborhood_human_transferred_maxSTRINGsimilarity.hdf5
TEXTMIN_HDF5_FILE=./regions_25kb_CT20220706.textmining_human_transferred_maxSTRINGsimilarity.hdf5
##################
GENOMIC_LOCATION=$1
##################

echo $GENOMIC_LOCATION

#find similar regions
echo "Finding Similar Regions"
time python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_combine_scores_CT20220707.py $SCALE_GZ_FILENAME $GENOMIC_LOCATION $DHS_COSINE_HDF5_FILE \
$TF_COSINE_HDF5_FILE $EPILOGOS_COSINE_HDF5_FILE $COEX_HDF5_FILE $COOCCUR_HDF5_FILE $DATABASE_HDF5_FILE $EXPER_HDF5_FILE $FUSION_HDF5_FILE $HOMOLOGY_HDF5_FILE $NEIGH_HDF5_FILE \
$TEXTMIN_HDF5_FILE $SORTED_SCALE_FILENAME $DHS_OOE_HDF5_FILE $TF_OOE_HDF5_FILE $EPILOGOS_OOE_HDF5_FILE $DHS_MAPPED_FILE $TF_MAPPED_FILE $EPILOGOS_MAPPED_FILE

echo "Done"
