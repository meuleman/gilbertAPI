import numpy as np
import os
from datetime import date
import sys
import re
import h5py

def ooe_calc(score_data, promoter_regions):
    #seperating score labels and score values
    score_data_mean_vals = score_data[:, 3::2].astype(float)
    score_data_max_vals = score_data[:, 4::2].astype(float)

    #finding the average mean and max of scores for each label
    expected_mean = np.mean(score_data_mean_vals[promoter_regions], axis=0)
    expected_max = np.mean(score_data_max_vals[promoter_regions], axis=0)

    #finding ooe for score mean and max
    ooe_mean = score_data_mean_vals / expected_mean
    ooe_max = score_data_max_vals / expected_max

    #combining sum and max
    ooe = np.hstack((ooe_mean, ooe_max))
    
    #set non-promoter regions to 0
    ooe_promoter_only = np.zeros(ooe.shape)
    ooe_promoter_only[promoter_regions] = ooe[promoter_regions]

    return ooe_promoter_only

def main(mapped_scores_filename, promoter_mapped_filename, out_path):
    #load score data
    score_data = []
    f = open(mapped_scores_filename, 'r')
    next(f)
    for line in f:
        score_data.append(re.split('\t|\|', line.strip()))
    f.close()
    score_data = np.array(score_data)

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

    #enrichment for scores in each region
    ooe = ooe_calc(score_data, promoter_regions)
    
    #save ooe array
    ooe_f = h5py.File(out_path, 'w')
    dset = ooe_f.create_dataset('ooe_arr', data=ooe)
    ooe_f.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
