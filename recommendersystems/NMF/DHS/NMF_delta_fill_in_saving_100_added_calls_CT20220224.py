import pandas as pd
import numpy as np
import gzip
import pickle
import os
import sys
from datetime import date
import random
from sklearn.decomposition import NMF, non_negative_factorization
import time

def saving_100_added_calls():
    #saving only the added delta values and biosample/DHS inds
    for i in range(2, 11):
        print(i, flush=True)
        #load npz file
        filename = '/home/ctrader/Hold/NMF/DHS/fill_in_multiple_at_a_time_saves/delta_fill_in_multiple_at_a_time_diffflat_sorteddiffbiosample_sorteddiffdhs_fillsize100_cycle'+str(i)+'_CT20220224.npz'
        cycle_np_file = np.load(filename)
    
        #load biosample and DHS inds
        cycle_sorted_diff_DHS_inds = cycle_np_file['sorted_diff_DHS_inds']
        cycle_sorted_diff_biosample_inds = cycle_np_file['sorted_diff_biosample_inds']
    
        #find number of DHSs
        n_DHSs = int(len(cycle_sorted_diff_DHS_inds)/733)
    
        #recreate flat indices
        cycle_sorted_diff_flat_inds = np.int64(cycle_sorted_diff_DHS_inds) + (n_DHSs * np.int64(cycle_sorted_diff_biosample_inds))
    
        #take the 100 added flat indices
        added_100_flat_inds = cycle_sorted_diff_flat_inds[-100:]
        del cycle_sorted_diff_flat_inds

        cycle_diff_flat_array = cycle_np_file['diff_flat_array']

        added_100_delta_vals = np.flip(cycle_diff_flat_array[added_100_flat_inds])
        added_100_DHS_inds = np.flip(cycle_sorted_diff_DHS_inds[-100:])
        added_100_biosample_inds = np.flip(cycle_sorted_diff_biosample_inds[-100:])

        outfile = '/home/ctrader/Hold/NMF/DHS/fill_in_multiple_at_a_time_saves/delta_fill_in_multiple_at_a_time_diffflat_sorteddiffbiosample_sorteddiffdhs_fillsize100_cycle'+str(i)+'_addedonly_CT20220224.npz'
        np.savez(outfile, delta_vals=added_100_delta_vals, DHS_inds=added_100_DHS_inds, biosample_inds=added_100_biosample_inds)
    
        del added_100_delta_vals
        del added_100_DHS_inds
        del added_100_biosample_inds
        del cycle_diff_flat_array
        del cycle_sorted_diff_DHS_inds
        del cycle_sorted_diff_biosample_inds

def main():
    saving_100_added_calls()

if __name__ == "__main__":
    main()
