import numpy as np
import os
from datetime import date
import sys
import h5py
import re

def ooe_calc(count_data, promoter_regions):
    #seperating count labels and count values
    count_data_vals = count_data[1:, 3:].astype(int)

    #finding the total number of counts for each label in promoter regions
    total_counts = np.sum(count_data_vals[promoter_regions], axis=0)

    #finding observed over expected for each label
    expected = total_counts / np.sum(total_counts)
    observed = (count_data_vals.T / np.sum(count_data_vals, axis=1)).T
    observed = np.where(np.isnan(observed), 0, observed)
    ooe = observed / expected

    #set non-promoter regions to 0
    ooe_promoter_only = np.zeros(ooe.shape)
    ooe_promoter_only[promoter_regions] = ooe[promoter_regions]
    return ooe_promoter_only

def main(mapped_counts_filename, promoter_mapped_filename, out_path):
    #load count data
    count_data = []
    f = open(mapped_counts_filename, 'r')
    for line in f:
        count_data.append(line.strip().split('\t'))
    f.close()
    count_data = np.array(count_data)

    #finding regions with promoters
    promoter_regions = []
    f = open(promoter_mapped_filename, 'r')
    region = 0
    for line in f:
        line = re.split('\t|;|\|', line.strip())
        if line[3] != '':
            promoter_regions.append(region)
        region += 1
    f.close()
    promoter_regions = np.array(promoter_regions)
    
    #enrichment for counts in each region
    ooe = ooe_calc(count_data, promoter_regions)
    
    #save ooe array
    ooe_f = h5py.File(out_path, 'w')
    dset = ooe_f.create_dataset('ooe_arr', data=ooe)
    ooe_f.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
