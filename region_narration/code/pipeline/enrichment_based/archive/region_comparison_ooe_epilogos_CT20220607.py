import numpy as np
import os
from datetime import date
import sys
import re
import pandas as pd

def ooe_calc(score_data, bed_files_path):
    # seperating score labels and score values
    score_val_inds = np.arange(3, len(score_data[0, :]), 3)
    score_data_labels = score_data[0, score_val_inds]
    score_data_sum_vals = score_data[:, score_val_inds+1].astype(float)
    score_data_max_vals = score_data[:, score_val_inds+2].astype(float)

    total_score_sum = np.array([], dtype=int)
    total_score_max = np.array([], dtype=int)
    for label in score_data_labels:
        label_path = os.path.join(bed_files_path, label+'.bed')
        label_data = np.array(pd.read_csv(label_path, sep='\t', usecols=[4], header=None, squeeze=True))
        total_score_sum = np.append(total_score_sum, np.sum(label_data))
        total_score_max = np.append(total_score_max, np.max(label_data))

    #finding ooe for score sums
    expected_sum = total_score_sum / np.sum(total_score_sum)
    observed_sum = (score_data_sum_vals.T / np.sum(score_data_sum_vals, axis=1)).T
    observed_sum = np.where(np.isnan(observed_sum), 0, observed_sum)
    ooe_sum = observed_sum / expected_sum

    #finding ooe for score max
    expected_max = total_score_max / np.max(total_score_max)
    observed_max = (score_data_max_vals.T / np.max(score_data_max_vals, axis=1)).T
    observed_max = np.where(np.isnan(observed_max), 0, observed_max)
    ooe_max = observed_max / expected_max

    #combining sum and max
    ooe = np.hstack((ooe_sum, ooe_max))

    return ooe

def main(mapped_scores_filename, bed_files_path, out_path):
    #load score data
    score_data = []
    f = open(mapped_scores_filename, 'r')
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
