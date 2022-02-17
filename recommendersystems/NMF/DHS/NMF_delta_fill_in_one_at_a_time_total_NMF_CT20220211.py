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

def fill_in_delta_one_at_a_time_total_NMF(decompBasis, decompMixture, dhs_arr, cycles):
    today=date.today()
    
    updated_dhs_array = np.copy(dhs_arr)
    Basis = np.copy(decompBasis)
    Mixture = np.copy(decompMixture)
    
    n_biosamples, n_biosample_e = Basis.shape
    n_DHS_e, n_DHSs = Mixture.shape
    
    changed_values = [] # biosample, DHS, changed to what value
    for i in range(cycles):
        t = time.time()
        print('Starting Cycle %s / %s' %(i+1, cycles), flush=True)
        #Generate scores for each interaction
        scores = np.dot(Basis, Mixture)
        #finding delta values between most recent presence/absence matrix and the reconstructed matrix        
        diff = scores - updated_dhs_array
        
        max_diff_abs = np.absolute(np.max(diff))
        min_diff_abs = np.absolute(np.min(diff))
        max_or_min = np.argmax([max_diff_abs, min_diff_abs])
#         print(time.time()-t)
        
        if max_or_min == 0:
            max_diff_flat_ind = np.argmax(diff)
            DHS_ind = max_diff_flat_ind%n_DHSs
            biosample_ind = max_diff_flat_ind//n_DHSs
            change_num = 1
            delta_val = np.max(diff)
        else:
            min_diff_flat_ind = np.argmin(diff)
            DHS_ind = min_diff_flat_ind%n_DHSs
            biosample_ind = min_diff_flat_ind//n_DHSs
            change_num = 0
            delta_val = np.min(diff)
#         print(time.time()-t)
        
        updated_dhs_array[biosample_ind, DHS_ind] = change_num
#         print(time.time()-t)
        #regenerate NMF decompositions
#         print()
#         print('NMF Decomposition:')
        
        changed_values.append([biosample_ind, DHS_ind, change_num, delta_val])

        new_Basis, new_Mixture, n_iter = non_negative_factorization(updated_dhs_array, W=Basis, H=Mixture, n_components=16, init='nndsvd', random_state=3, max_iter=400, verbose=1)
        
        Basis = new_Basis
        Mixture = new_Mixture
        
        np.savez('/home/ctrader/Hold/NMF/DHS/fill_in_one_at_a_time_saves/one_at_a_time_cycle_'+str(i)+'_basis_mixture_total_NMF_CT'+today.strftime("%Y%m%d")+'.npy', basis=Basis, mixture=Mixture)
        print(time.time()-t)
    changed_values = np.array(changed_values)
    np.save('/home/ctrader/Hold/NMF/DHS/fill_in_one_at_a_time_saves/one_at_a_time_changed_interactions_total_NMF_'+str(cycles)+'_cycles_CT'+today.strftime("%Y%m%d")+'.npy', changed_values)
    return updated_dhs_array, Basis, Mixture, changed_values

def main():
    #read in 733 sample x 3.5+ million DHS presence/absence matrix
    dhs = pd.read_table('/home/ctrader/Hold/DHS/dat_bin_FDR01_hg38.txt.gz', header=None).T
    #Loading biosample embeddings
    fn = gzip.GzipFile('/home/ctrader/Hold/DHS/2018-06-08NC16_NNDSVD_Basis.npy.gz', "r")
    decomp_Basis = np.load(fn)
    #Loading DHS embeddings
    fn = gzip.GzipFile('/home/ctrader/Hold/DHS/2018-06-08NC16_NNDSVD_Mixture.npy.gz', "r")
    decomp_Mixture = np.load(fn)
    dhs_array = dhs.to_numpy()
    updated_dhs_array, Basis, Mixture, changed_values = fill_in_delta_one_at_a_time_total_NMF(decompBasis=decomp_Basis, decompMixture=decomp_Mixture, dhs_arr=dhs_array, cycles=1000)

if __name__ == "__main__":
    main()
