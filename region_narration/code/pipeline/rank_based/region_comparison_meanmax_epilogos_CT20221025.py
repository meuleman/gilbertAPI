import numpy as np
import os
from datetime import date
import sys
import re
import h5py

def reorder_enrichments(score_data):
    #seperating score labels and score values
    score_data_mean_vals = score_data[:, 3::2].astype(float)
    score_data_max_vals = score_data[:, 4::2].astype(float)

    #combining sum and max
    reordered = np.hstack((score_data_mean_vals, score_data_max_vals))

    return reordered

def main(mapped_scores_filename, out_path):
    #load score data
    score_data = []
    f = open(mapped_scores_filename, 'r')

    next(f)
    for line in f:
        score_data.append(re.split('\t|\|', line.strip()))
    f.close()
    score_data = np.array(score_data)
    
    #enrichment for scores in each region
    reordered = reorder_enrichments(score_data)
    
    #save ooe array
    ooe_f = h5py.File(out_path, 'w')
    dset = ooe_f.create_dataset('base_arr', data=reordered)
    ooe_f.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
