#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import h5py
import csv
import numpy as np

# Open the CSV file
gene_length_list = []

# Open the CSV file
with open('/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/data_results/edited_gene_name_length_idx_data.csv', 'r') as f:
    # Skip the header row
    next(f)

     # Read the remaining rows
    reader = csv.reader(f)
    for row in reader:
        gene_length_list.append(row[3])

gene_length_list = [float(x) for x in gene_length_list]
print("Finished reading gene length information.")

# Open the h5py file
with h5py.File('/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files/makefiles/192k_biosamples_57k_genes_extractedMatrix_LW20230125.h5', 'r+') as hf:
    # Get the matrix from the file
    matrix = hf['matrix']
    matrix = np.array(matrix)
    matrix = matrix.astype(np.float32)
    print("Finished reading matrix.")

    # Divide all the values in the matrix by the corresponding values in the gene_length_list
    for i, length in enumerate(gene_length_list):
        if length == 0:
            continue
        matrix[i, :] = matrix[i, :] / length
    print("Finished normalizing matrix.")
     

    # Save the new matrix in the h5py file
    hf.create_dataset('corrected_matrix', data=matrix)

    # Create a dataset for the gene_length_list
    hf.create_dataset('gene_length', data=gene_length_list)

    # Delete the original matrix from the file
    #del hf['matrix']

    print("Finished saving corrected matrix and gene length information.")
