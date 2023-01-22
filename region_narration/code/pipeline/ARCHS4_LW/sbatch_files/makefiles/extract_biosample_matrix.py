#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import json
import pickle
import h5py

# Open ARCHS4 v2.1 Data
infile = "/home/ajoshi/archs4/data/archs4_gene_human_v2.1.h5"

# Open the JSON file
# root = "/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/results/"
# with open(root + 'allSelected_uniqueSampleIndices_list.json') as json_file:
#     # Load the JSON data from the file
#     sampleIdxs_list = json.load(json_file)

# Read ARCHS4 Matrix
h5py_file = h5py.File(infile, 'r')

# biosample_metadata = h5py_file[:, sampleIdxs_list]
# gene_metadata = h5py_file[geneIdxs_list, :]
biosample_metadata = h5py_file['biosample_metadata_path'][cols]
gene_metadata = h5py_file['gene_metadata_path'][rows]

biosample_gene_metadata = h5py_file[geneIdxs_list, sampleIdxs_list]
h5py_file.close()

# Output 
out_f = h5py.File(root + "biosample_gene_metadata_extractedMatrix.hdf5", 'w')
out_f.create_dataset('matrix', data=biosample_gene_metadata)

# out_f.create_dataset('biosample_metadata_extractedMatrix', data=biosample_metadata)
# out_f.create_dataset('gene_metadata_extractedMatrix', data=gene_metadata)
out_f.create_dataset('biosample_metadata', data=biosample_metadata)
out_f.create_dataset('gene_metadata', data=gene_metadata)

out_f.close()
