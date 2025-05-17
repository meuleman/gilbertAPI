#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np 
import h5py
import csv
import pandas as pd 

infile = "/home/ctrader/region_similarity/data/archs4_human_matrix_v11.h5"
hdf5_file = h5py.File(infile, 'r')
data = hdf5_file['data']
expression = data['expression'][:]
human_expressed = np.divide(np.count_nonzero(expression > 0, axis=1),np.count_nonzero(~np.isnan(expression), axis=1))

meta = hdf5_file['meta']
genes = meta['genes']
symbols = genes['gene_symbol'][:]

symbols_decoded = (s.decode('utf-8', 'ignore') for s in symbols) # Convert from bytes to string 
dictionary = dict(zip(symbols_decoded, human_expressed))

hdf5_file.close()

print("dictionary", dictionary)
# Write gene name & averages to file
with open("gene_exp_percentage.tsv", 'w') as f: 
    for key, value in dictionary.items(): 
        f.write('%s\t%f\n' % (key, value))
f.close()       
