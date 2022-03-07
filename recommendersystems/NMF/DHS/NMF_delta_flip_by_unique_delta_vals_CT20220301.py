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

def flip_delta(Basis, Mixture, dhs_arr, cycles, flip_size):
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

        del diff

        #finding the unique delta values
        diff_flat_unique = np.unique(diff_flat, return_inverse=True)

        #indices to set to 0 - taking the bottom flip_size unique delta values and flipping corresponding calls
        #to_0_inds = sorted_diff_flat[:flip_size]
        to_0_inds = np.where(diff_flat_unique[1] < flip_size)[0]
        dhs_arr_flat[to_0_inds] = 0

        #indices to set to 1 - taking the top flip_size unique delta values and flipping corresponding calls
        #to_1_inds = sorted_diff_flat[-flip_size:]
        to_1_inds = np.where(diff_flat_unique[1] > np.max(diff_flat_unique[1]) - flip_size)[0]
        dhs_arr_flat[to_1_inds] = 1

        del diff_flat_unique

        print('Sorting...', flush=True)
        sorted_diff_flat = np.argsort(diff_flat)

        dhs_arr = dhs_arr_flat.reshape((n_biosamples, n_DHSs))
        
        #regenerate NMF decompositions
        print('NMF Decomposition:', flush=True)
        Basis, Mixture, n_iter = non_negative_factorization(dhs_arr, W=np.float64(Basis), H=np.float64(Mixture), n_components=16, init='custom', max_iter=400, verbose=1)

        Basis = np.float32(Basis)
        Mixture = np.float32(Mixture)
        
        outfile = '/home/ctrader/Hold/NMF/DHS/fill_in_multiple_at_a_time_saves/delta_fill_in_multiple_at_a_time_booldhsarr_basis_mixture_fillsize'+str(flip_size)+'unique_cycle'+str(i+1)+'_CT'+today.strftime("%Y%m%d")+'.npz'
        np.savez(outfile, new_dhs_array=dhs_arr.astype(bool), new_decomp_basis=Basis, new_decomp_mixture=Mixture)

        diff_flat = np.float32(diff_flat)
        
        sorted_diff_flat_biosample_inds = np.int16(sorted_diff_flat // n_DHSs)
        sorted_diff_flat_DHS_inds = np.int32(sorted_diff_flat % n_DHSs)

        outfile = '/home/ctrader/Hold/NMF/DHS/fill_in_multiple_at_a_time_saves/delta_fill_in_multiple_at_a_time_diffflat_sorteddiffbiosample_sorteddiffdhs_fillsize'+str(flip_size)+'unique_cycle'+str(i+1)+'_CT'+today.strftime("%Y%m%d")+'.npz'
        np.savez(outfile, diff_flat_array=diff_flat, sorted_diff_DHS_inds=sorted_diff_flat_DHS_inds, sorted_diff_biosample_inds=sorted_diff_flat_biosample_inds, inds_to_1=to_1_inds, inds_to_0=to_0_inds)

        del sorted_diff_flat
        del diff_flat
        del sorted_diff_flat_DHS_inds
        del sorted_diff_flat_biosample_inds
        del to_1_inds
        del to_0_inds

    return    

def main(argv1, argv2):
    #read in 733 sample x 3.5+ million DHS presence/absence matrix
    dhs = pd.read_table('/home/ctrader/Hold/DHS/dat_bin_FDR01_hg38.txt.gz', header=None).T
    
    #Loading biosample embeddings
    fn = gzip.GzipFile('/home/ctrader/Hold/DHS/2018-06-08NC16_NNDSVD_Basis.npy.gz', "r")
    decomp_Basis = np.float32(np.load(fn))
    
    #Loading DHS embeddings
    fn = gzip.GzipFile('/home/ctrader/Hold/DHS/2018-06-08NC16_NNDSVD_Mixture.npy.gz', "r")
    decomp_Mixture = np.float32(np.load(fn))
    
    dhs_array = dhs.to_numpy()
    del dhs
    
    cycles = int(argv1)
    flip_size = int(argv2)
    
    flip_delta(decomp_Basis, decomp_Mixture, dhs_array, cycles, flip_size)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
