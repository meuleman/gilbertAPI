#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import csv
import json
import pickle
import h5py
import numpy as np

# Open ARCHS4 v2.1 Data
infile = "/home/ajoshi/archs4/data/archs4_gene_human_v2.1.h5"

# Open the JSON file - Selected Samples
root = "/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/results/"
with open(root + 'allSelected_uniqueSampleIdx_list_LW20230121.json') as json_file:
    # Load the JSON data from the file
    sampleIdxs_list = json.load(json_file)

# Open the csv file - Selected Genes 
gene_names = []
gene_idxs = []

with open("/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/data_results/edited_gene_name_length_idx_data.csv", 'r') as file:
    reader = csv.reader(file)
    # Skip the header row
    next(reader)
    for row in reader:
        gene_names.append(str(row[1]))
        gene_idxs.append(int(row[4]))

# Get the selected samples - names 
pickle_in = open("/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/Analysis/afterReadFilteration_CategorySampleSizes/archs4_sampleNames_list.pkl", 'rb')
sample_names_list = pickle.load(pickle_in)
selected_sample_names = [sample_names_list[i] for i in sampleIdxs_list]
pickle_in.close()

# Read ARCHS4 Matrix
with h5py.File(infile, 'r') as h5py_file:
    data = h5py_file['data']
    matrix = data['expression'][:]

# Extract the rows and columns specified by the indexes
extracted_matrix = matrix[gene_idxs, :][:, sampleIdxs_list]

# Create a new H5 file
with h5py.File('extracted_matrix.h5', 'w') as f:
    
    # Create a dataset for the matrix
    dset = f.create_dataset('matrix', data=extracted_matrix)
     # Create a dataset for the sample names
    dset2 = f.create_dataset('sample_names', (len(selected_sample_names),), dtype=h5py.special_dtype(vlen=str))
    dset2[:] = selected_sample_names
    # Create a dataset for the gene names
    dset3 = f.create_dataset('gene_names', (len(gene_names),), dtype=h5py.special_dtype(vlen=str))
    dset3[:] = gene_names
