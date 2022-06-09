import numpy as np
import os
from datetime import date
import sys
import re

def ooe_calc(score_data, bed_files_path):
    #seperating score labels and score values
    score_data_mean_vals = score_data[:, 3::2].astype(float)
    score_data_max_vals = score_data[:, 4::2].astype(float)

    #finding the average mean and max of scores for each label
    expected_mean = np.mean(score_data_mean_vals, axis=0)
    expected_max = np.mean(score_data_max_vals, axis=0)

    #finding ooe for score mean and max
    ooe_mean = score_data_mean_vals / expected_mean
    ooe_max = score_data_max_vals / expected_max

    #combining sum and max
    ooe = np.hstack((ooe_mean, ooe_max))

    return ooe

def main(mapped_scores_filename, bed_files_path, out_path):
    #load score data
    score_data = []
    f = open(mapped_scores_filename, 'r')

    next(f)
    for line in f:
        score_data.append(re.split('\t|\|', line.strip()))
    f.close()
    score_data = np.array(score_data)
    
    #enrichment for scores in each region
    ooe = ooe_calc(score_data, bed_files_path)
    
    #save ooe array
    np.save(out_path, ooe)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
