import numpy as np
import os
from datetime import date
import sys

def ooe_calc(count_data, bed_files_path):
    #seperating count labels and count values
    count_val_inds = np.arange(4, len(count_data[0, :]), 2)
    count_data_vals = count_data[:, count_val_inds].astype(int)
    count_data_labels = count_data[0, count_val_inds-1]

    #finding the total number of counts for each label
    total_counts = np.array([], dtype=int)
    for label in count_data_labels:
        label_path = os.path.join(bed_files_path, label+'.bed')
        f = open(label_path, 'r')
        total_counts = np.append(total_counts, len(f.readlines()))

    #finding observed over expected for each label
    #expected = total_counts / len(count_data_vals)
    #ooe = count_data_vals/expected
    expected = total_counts / np.sum(total_counts)
    observed = (count_data_vals.T / np.sum(count_data_vals, axis=1)).T
    observed = np.where(np.isnan(observed), 0, observed)
    ooe = observed / expected
    return ooe

def main(mapped_counts_filename, bed_files_path, out_path):
    #load count data
    count_data = []
    f = open(mapped_counts_filename, 'r')
    for line in f:
        count_data.append(line.strip().split('\t'))
    f.close()
    count_data = np.array(count_data)
    
    #enrichment for counts in each region
    ooe = ooe_calc(count_data, bed_files_path)
    
    #save ooe array
    np.save(out_path, ooe)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
