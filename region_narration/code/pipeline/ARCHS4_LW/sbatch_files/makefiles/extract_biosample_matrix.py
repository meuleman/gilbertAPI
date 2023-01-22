#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import json
import pickle
import h5py

# Open ARCHS4 v2.1 Data
infile = "/home/ajoshi/archs4/data/archs4_gene_human_v2.1.h5"

# Open the JSON file
root = "/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/results/"

import csv

with open("/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/data_results/gene_name_length_idx_data.csv", 'r') as file:
    reader = csv.reader(file)
    # Skip the header row
    next(reader)
    gene_names = [str(row[1]) for row in reader]
    gene_idxs = [int(row[-1]) for row in reader]


with open(root + 'allSelected_uniqueSampleIndices_list.json') as json_file:
    # Load the JSON data from the file
    sampleIdxs_list = json.load(json_file)

pickle_in = open("/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/Analysis/afterReadFilteration_CategorySampleSizes/archs4_sampleNames_list.pkl", 'rb')
gene_names_list = pickle.load(pickle_in)
list3 = [gene_names_list[i] for i in sampleIdxs_list]

# Read ARCHS4 Matrix
with h5py.File(infile, 'r') as h5py_file:
    data = h5py_file['data']
    expression = data['expression'][:]

biosample_metadata = h5py_file[:, sampleIdxs_list]

biosample_metadata = h5py_file['biosample_metadata_path'][list3]

extractedMatrix = expression[gene_idxs, sampleIdxs_list]
gene_metadata = expression['gene_metadata_path'][gene_names]
h5py_file.close()

# Output 
out_f = h5py.File(root + "gene_extractedMatrix.hdf5", 'w')
out_f.create_dataset('matrix', data=extractedMatrix)

# out_f.create_dataset('biosample_metadata_extractedMatrix', data=biosample_metadata)
# out_f.create_dataset('gene_metadata_extractedMatrix', data=gene_metadata)
#out_f.create_dataset('biosample_metadata', data=biosample_metadata)
out_f.create_dataset('gene_metadata', data=gene_metadata)

out_f.close()
