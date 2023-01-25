#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import h5py
import csv
import qnorm 

# Open the h5py file
with h5py.File('/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/makefiles/192k_biosamples_57k_genes_correctedMatrix_LW20230123.h5', 'r+') as hf:
    # Get the matrix from the file
    matrix = hf['new_matrix']

    # Quantile normalize 
    qnorm_matrix = qnorm.quantile_normalize(np.log2(1+np.array(matrix)))
        matrix[i, :] = matrix[i, :] / length

    # Save the new matrix in the h5py file
    hf.create_dataset('matrix', data=qnorm_matrix)

    # Delete the original matrix from the file
    del hf['new_matrix']
