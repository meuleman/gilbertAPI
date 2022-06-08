import numpy as np
import sys
import os

def combine_scores():
    

def main(region_of_interest_path, DHS_npz_path, TF_npz_path):
    #find the region of interest
    f = open(region_of_interest_path, 'r')
    roi = int(f.readline())

    #####Load rank and ind arrays (np.load('')['this array'][roi]
    np.load(DHS_npz_path)['

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
