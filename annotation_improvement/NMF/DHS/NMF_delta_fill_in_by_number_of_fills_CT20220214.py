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

def fill_in_delta(Basis, Mixture, dhs_arr, cycles, flip_size):
    today=date.today()
    
    n_biosamples, n_biosample_e = Basis.shape
    n_DHS_e, n_DHSs = Mixture.shape

    for i in range(cycles):
        print('Starting Cycle %s / %s' %(i+1, cycles), flush=True)

        #Generate scores for each call
        print('Reconstructing DHS Matrix...', flush=True)
        scores = np.dot(Basis, Mixture)

        #finding delta values between most recent presence/absence matrix and the reconstructed matrix
        diff = scores - dhs_arr 

        print('Updating presence/absence matrix...', flush=True)
        dhs_arr_flat = dhs_arr.flatten()
        diff_flat = diff.flatten()

        print('Sorting...', flush=True)
        sorted_diff_flat = np.argsort(diff_flat)

        #indices to set to 0
        to_0_inds = sorted_diff_flat[:flip_size]
        dhs_arr_flat[to_0_inds] = 0

        #indices to set to 1
        to_1_inds = sorted_diff_flat[-flip_size:]
        dhs_arr_flat[to_1_inds] = 1

        dhs_arr = dhs_arr_flat.reshape((n_biosamples, n_DHSs))
        
        #regenerate NMF decompositions
        print('NMF Decomposition:', flush=True)
        Basis, Mixture, n_iter = non_negative_factorization(dhs_arr, W=Basis, H=Mixture, n_components=16, init='custom', max_iter=400, verbose=1)
        
        outfile = '/home/ctrader/Hold/NMF/DHS/fill_in_multiple_at_a_time_saves/delta_fill_in_multiple_at_a_time_dhsarr_basis_mixture_fillsize'+str(flip_size)+'_cycle'+str(i+1)+'_CT'+today.strftime("%Y%m%d")+'.npz'
        np.savez(outfile, new_dhs_array=dhs_arr, new_decomp_basis=Basis, new_decomp_mixture=Mixture)

        outfile = '/home/ctrader/Hold/NMF/DHS/fill_in_multiple_at_a_time_saves/delta_fill_in_multiple_at_a_time_diffflat_sorteddiffflat_fillsize'+str(flip_size)+'_cycle'+str(i+1)+'_CT'+today.strftime("%Y%m%d")+'.npz'
        np.savez(outfile, diff_flat_array=diff_flat, sorted_diff_flat_array=sorted_diff_flat)
    return dhs_arr, Basis, Mixture, diff_flat, sorted_diff_flat    

def main(argv1, argv2):
    #read in 733 sample x 3.5+ million DHS presence/absence matrix
    dhs = pd.read_table('/home/ctrader/Hold/DHS/dat_bin_FDR01_hg38.txt.gz', header=None).T
    #Loading biosample embeddings
    fn = gzip.GzipFile('/home/ctrader/Hold/DHS/2018-06-08NC16_NNDSVD_Basis.npy.gz', "r")
    decomp_Basis = np.load(fn)
    #Loading DHS embeddings
    fn = gzip.GzipFile('/home/ctrader/Hold/DHS/2018-06-08NC16_NNDSVD_Mixture.npy.gz', "r")
    decomp_Mixture = np.load(fn)
    dhs_array = dhs.to_numpy()
    cycles = int(argv1)
    flip_size = int(argv2)
    updated_dhs_array, Basis, Mixture, changed_values = fill_in_delta(decomp_Basis, decomp_Mixture, dhs_array, cycles, flip_size)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
