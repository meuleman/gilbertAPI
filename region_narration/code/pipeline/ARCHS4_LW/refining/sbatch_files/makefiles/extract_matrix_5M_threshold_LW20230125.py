#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import csv
import pickle
import h5py
import numpy as np

# Open the samples idxs above 5M
root = "/home/lwu/ARCHS4/Tests/"
pickle_in = open(root + '5M_threshold_sample_idxs_list.pkl', 'rb')
sampleIdxs_5M_list = pickle.load(pickle_in)
pickle_in.close()

# Open the h5py file
with h5py.File('/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/makefiles/192k_biosamples_57k_genes_correctedMatrix_LW20230125.h5', 'r+') as hf:
    # Get the matrix from the file
    matrix = hf['corrected_matrix'][:]

    extracted_matrix = matrix[:][:, sampleIdxs_5M_list]
    
    # Save the new matrix in the h5py file
    hf.create_dataset('threshold_matrix', data=extracted_matrix)

