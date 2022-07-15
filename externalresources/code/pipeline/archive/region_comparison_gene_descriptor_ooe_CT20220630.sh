#!/bin/bash
#SBATCH --job-name=region_comparison
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10000
#SBATCH --error=/home/ctrader/region_similarity/slurm_out_error/gene_descriptor_ooe.err
#SBATCH --output=/home/ctrader/region_similarity/slurm_out_error/gene_descriptor_ooe.out
pwd; hostname; date

module load anaconda
module load bedops

eval "$(conda shell.bash hook)"
conda activate lightfm

cd /home/ctrader/region_similarity/data

#####Variables####
EPILOGOS_COSINE_HDF5_FILE=./regions_25kb_CT20220623.epilogoscosine.hdf5
TF_COSINE_HDF5_FILE=./regions_25kb_CT20220623.TFcosine.hdf5
DHS_COSINE_HDF5_FILE=./regions_25kb_CT20220623.DHScosine.hdf5

GENES_BED=./mart_export_protein_coding_genes_unique_cleaned_CT20220628.bed
GENES_MAPPED_BED=./regions_25kb_CT20220628.cleanedgenesmapped.bed
GENE_DESCRIPTOR_OOE=./regions_25kb_CT20220630.genedescriptorooe.hdf5
GENE_DESCRIPTORS=./gene_descriptors_CT20220630.npy

NUM_CPUS=10
##################

#Calculating OOE
echo "Calculating OOE for Gene Descriptors"
python /home/ctrader/AG/AugmentedGenomics/externalresources/code/pipeline/region_comparison_gene_descriptor_ooe_CT20220630.py $DHS_COSINE_HDF5_FILE $TF_COSINE_HDF5_FILE $EPILOGOS_COSINE_HDF5_FILE $GENES_BED $GENES_MAPPED_BED $GENE_DESCRIPTOR_OOE $GENE_DESCRIPTORS $NUM_CPUS

echo "Done"
