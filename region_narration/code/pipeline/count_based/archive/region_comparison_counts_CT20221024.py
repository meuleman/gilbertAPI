import numpy as np
import os
from datetime import date
import sys
import h5py

def main(mapped_counts_filename, out_path):
    #load count data
    count_data = []
    f = open(mapped_counts_filename, 'r')
    for line in f:
        count_data.append(line.strip().split('\t'))
    f.close()
    count_data = np.array(count_data)
    
    #seperating count labels and count values
    count_data_vals = count_data[1:, 3:].astype(int)

    #save count array
    count_f = h5py.File(out_path, 'w')
    dset = count_f.create_dataset('count_arr', data=count_data_vals)
    count_f.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
