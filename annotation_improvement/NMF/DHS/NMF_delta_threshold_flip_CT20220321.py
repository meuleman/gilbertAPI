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

def flip_delta(Basis, Mixture, dhs_arr, cycles, to0_thresh, to1_thresh):
    today=date.today()
    
    n_biosamples, n_biosample_e = Basis.shape
    n_DHS_e, n_DHSs = Mixture.shape

    for i in range(cycles):
        print('Starting Cycle %s / %s' %(i+1, cycles), flush=True)

        #Generate scores for each call
        print('Reconstructing DHS Matrix...', flush=True)
        scores = np.dot(Basis, Mixture)

        #finding delta values between most recent presence/absence matrix and the reconstructed matrix
        diff = np.float32(scores - dhs_arr)

        del scores

        print('Updating presence/absence matrix...', flush=True)
        diff_flat = np.float32(diff.flatten())
        
        #flipping calls with delta values smaller than to0 threshold to 0
        to_0_inds = np.where(diff_flat <= to0_thresh)[0]
        to_0_biosample_inds = np.int16(to_0_inds // n_DHSs)
        to_0_DHS_inds = np.int32(to_0_inds % n_DHSs)
        del to_0_inds
        dhs_arr = np.where(diff <= to0_thresh, 0, dhs_arr)

        #flipping calls with delta values larger than to1 threshold to 1
        to_1_inds = np.where(diff_flat >= to1_thresh)[0]
        to_1_biosample_inds = np.int16(to_1_inds // n_DHSs)
       	to_1_DHS_inds =	np.int32(to_1_inds % n_DHSs)
        del to_1_inds
        dhs_arr	= np.where(diff >= to1_thresh, 1, dhs_arr)

        del diff

        print('Sorting...', flush=True)
        sorted_diff_flat = np.argsort(diff_flat)
        sorted_diff_flat_biosample_inds = np.int16(sorted_diff_flat // n_DHSs)
        sorted_diff_flat_DHS_inds = np.int32(sorted_diff_flat % n_DHSs)
        del sorted_diff_flat

        #regenerate NMF decompositions
        print('NMF Decomposition:', flush=True)
        Basis, Mixture, n_iter = non_negative_factorization(dhs_arr, W=Basis, H=Mixture, n_components=16, init='custom', max_iter=400, verbose=1)

        outfile = '/home/ctrader/Hold/NMF/DHS/fill_in_multiple_at_a_time_saves/delta_threshold_flip_booldhsarr_basis_mixture_to0thresh'+str(to0_thresh)+'_to1thresh'+str(to1_thresh)+'_cycle'+str(i+1)+'_CT'+today.strftime("%Y%m%d")+'.npz'
        np.savez(outfile, new_dhs_array=dhs_arr.astype(bool), new_decomp_basis=np.float32(Basis), new_decomp_mixture=np.float32(Mixture))
        
        outfile = '/home/ctrader/Hold/NMF/DHS/fill_in_multiple_at_a_time_saves/delta_threshold_flip_diffflat_sortedDHS_sortedbiosample_to0inds_to1inds_to0thresh'+str(to0_thresh)+'_to1thresh'+str(to1_thresh)+'_cycle'+str(i+1)+'_CT'+today.strftime("%Y%m%d")+'.npz'
        np.savez(outfile, diff_flat_array=diff_flat, sorted_diff_DHS_inds=sorted_diff_flat_DHS_inds, sorted_diff_biosample_inds=sorted_diff_flat_biosample_inds, biosample_inds_to_1=to_1_biosample_inds, DHS_inds_to_1=to_1_DHS_inds, biosample_inds_to_0=to_0_biosample_inds, DHS_inds_to_0=to_0_DHS_inds)

        del diff_flat
        del sorted_diff_flat_DHS_inds
        del sorted_diff_flat_biosample_inds
        del to_1_biosample_inds
        del to_1_DHS_inds
        del to_0_biosample_inds
        del to_0_DHS_inds

    return    

def main(argv1, argv2, argv3):
    #read in 733 sample x 3.5+ million DHS presence/absence matrix
    dhs = pd.read_table('/home/ctrader/Hold/DHS/dat_bin_FDR01_hg38.txt.gz', header=None).T
    
    #Loading biosample embeddings
    fn = gzip.GzipFile('/home/ctrader/Hold/DHS/2018-06-08NC16_NNDSVD_Basis.npy.gz', "r")
    decomp_Basis = np.load(fn)
    
    #Loading DHS embeddings
    fn = gzip.GzipFile('/home/ctrader/Hold/DHS/2018-06-08NC16_NNDSVD_Mixture.npy.gz', "r")
    decomp_Mixture = np.load(fn)
    
    dhs_array = dhs.to_numpy()
    del dhs
    
    cycles = int(argv1)
    to0_thresh = float(argv2)
    to1_thresh = float(argv3)
    
    flip_delta(decomp_Basis, decomp_Mixture, dhs_array, cycles, to0_thresh, to1_thresh)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
