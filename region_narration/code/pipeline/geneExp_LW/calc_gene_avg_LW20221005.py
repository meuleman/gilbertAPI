#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np 
import h5py
import csv
import pandas as pd 

infile = "/home/ctrader/region_similarity/data/archs4_human_matrix_v11.h5"
hdf5_file = h5py.File(infile, 'r')
data = hdf5_file['data']
expression = data['expression'][:]
human_avgs = np.mean(expression, axis = 1)

meta = hdf5_file['meta']
genes = meta['genes']
symbols = genes['gene_symbol'][:]

symbols_decoded = (s.decode('utf-8', 'ignore') for s in symbols) # Convert from bytes to string 
dictionary = dict(zip(symbols_decoded, human_avgs))

############################################ GMT DATASETS ####################################################### 

gene_set1_file = "/home/ctrader/region_similarity/data/h.all.v7.5.1.symbols.gmt"
gene_set2_file = "/home/ctrader/region_similarity/data/c1.all.v7.5.1.symbols.gmt"
gene_set3_file = "/home/ctrader/region_similarity/data/c2.all.v7.5.1.symbols.gmt"
row_avgs_list = []

# Lists for dataset additional details (ie. no data genes)
no_data_name = [] # genome names that has no-data genes 
notNaN_row_count = [] # genome count for non-NaN genes (inital dataset)
no_data_genes = [] # genes that have no data 
result_row_count = []  # genome count after filtering no-data genes

# Brute force in quickly gathering names (could revise)
with open(gene_set1_file, 'r') as file:
    csvreader = csv.reader(file, delimiter='\t')
    name = [row.split()[0] for row in file] 

file.close()
    
with open(gene_set1_file, 'r') as file:
    csvreader = csv.reader(file, delimiter='\t')

    for row in csvreader:

        row_sum = 0 # the sum of the averages in a row 
        row_avg = 0 # the avg of the averages in a row
        found_count = 0 # count of sucessfully found genes
        indv_no_data = [] # saves each row's no-data genes 
        
        notNaN_row_count.append(len(row)-2)
        
        for idx in range(2, len(row)):
            item = row[idx]

            if item in dictionary: 
                row_sum = row_sum + dictionary[item]
                found_count = found_count + 1
                
            else: 
                if row[0] not in no_data_name: 
                    no_data_name.append(row[0])
                indv_no_data.append(item)
        
        no_data_genes.append(indv_no_data)
        result_row_count.append(found_count)
        row_avg = row_sum/found_count
        row_avgs_list.append(row_avg) 
        
new_dictionary = dict(zip(name, row_avgs_list))

file.close()

# Write gene name & averages to file
with open("file1_avgs.tsv", 'w') as f: 
    for key, value in new_dictionary.items(): 
        f.write('%s\t%f\n' % (key, value))
f.close()       

# Write no-data details to file
with open('file1_details.tsv',"w") as f:
    for (no_data_name, notNaN_row_count, no_data_genes,result_row_count) in zip(no_data_name, notNaN_row_count, no_data_genes,result_row_count):
        f.write("{0}\t{1}\t{2}\t{3}\n".format(no_data_name, notNaN_row_count, no_data_genes,result_row_count))
f.close()   


# In[17]:


no_data_name = []
notNaN_row_count = []
no_data_genes = []
result_row_count = []    
with open(gene_set1_file, 'r') as file:
    csvreader = csv.reader(file, delimiter='\t')

    for row in csvreader:

        row_sum = 0 # the sum of the averages in a row 
        row_avg = 0 # the avg of the averages in a row
        found_count = 0
        indv_no_data = []
        
        print("row size:", len(row)-2)
        notNaN_row_count.append(len(row)-2)

        for idx in range(2, len(row)):
            item = row[idx]
            #print("item", item)

            if item in dictionary: 
                #print("row_sum", row_sum)
                row_sum = row_sum + dictionary[item]
                found_count = found_count + 1
                
        
            else: 
                if row[0] not in no_data_name: 
                    no_data_name.append(row[0])
                    
                
                print("name", row[0])
                
                print("false", item)
                
                indv_no_data.append(item)
        no_data_genes.append(indv_no_data)
        print("found_count", found_count)
        result_row_count.append(found_count)
        print()
        row_avg = row_sum/found_count
        row_avgs_list.append(row_avg) 
        #break

with open('c1.all.v7.5.1.symbols.tsv',"w") as f:
    for (no_data_name, notNaN_row_count, no_data_genes,result_row_count) in zip(no_data_name, notNaN_row_count, no_data_genes,result_row_count):
        f.write("{0}\t{1}\t{2}\t{3}\n".format(no_data_name, notNaN_row_count, no_data_genes,result_row_count))


# In[ ]:




