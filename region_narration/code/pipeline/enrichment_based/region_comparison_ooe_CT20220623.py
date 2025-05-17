import numpy as np
import os
from datetime import date
import sys
import h5py

def ooe_calc(count_data):
    #seperating count labels and count values
    count_data_vals = count_data[1:, 3:].astype(int)

    #finding the total number of counts for each label
    total_counts = np.sum(count_data_vals, axis=0)

    #finding observed over expected for each label
    #expected = total_counts / len(count_data_vals)
    #ooe = count_data_vals/expected
    expected = total_counts / np.sum(total_counts)
    observed = (count_data_vals.T / np.sum(count_data_vals, axis=1)).T
    observed = np.where(np.isnan(observed), 0, observed)
    ooe = observed / expected
    return ooe

def main(mapped_counts_filename, out_path):
    #load count data
    count_data = []
    f = open(mapped_counts_filename, 'r')
    for line in f:
        count_data.append(line.strip().split('\t'))
    f.close()
    count_data = np.array(count_data)
    
    #enrichment for counts in each region
    ooe = ooe_calc(count_data)
    
    #save ooe array
    ooe_f = h5py.File(out_path, 'w')
    dset = ooe_f.create_dataset('ooe_arr', data=ooe)
    ooe_f.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
