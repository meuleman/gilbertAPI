import numpy as np
from functools import partial
from multiprocessing import Pool
from contextlib import closing
import sys
import h5py
from collections import Counter
import re
from scipy.stats import fisher_exact

def find_significance_for_region(total_term_counts, region):
    #for a single region, check if each term that occurs is significant and return mask array
    region_counts = term_full_counts[region]
    terms_inds_to_test = np.where(region_counts!=0)[0]
    sig_pa = np.array([False]*len(region_counts))
    for term in terms_inds_to_test:
        term_count_in_region = region_counts[term]
        not_term_count_in_region = np.sum(region_counts) - term_count_in_region
        term_count_in_not_region = np.sum(term_full_counts[:, term]) - term_count_in_region
        not_term_count_in_not_region = total_term_counts - term_count_in_region - not_term_count_in_region - term_count_in_not_region
        #calculating p value from fisher exact test
        conting = [[term_count_in_region, not_term_count_in_region], 
               [term_count_in_not_region, not_term_count_in_not_region]]
        odds, pval = fisher_exact(conting, alternative='greater')
        #bonferroni correction
        if pval < (0.05 / len(terms_inds_to_test)):
            sig_pa[term] = True
    return sig_pa

def main(descriptorooe_hdf5_path, out_path, num_cpus):
    #load data
    global term_full_counts
    descriptor_hdf5 = h5py.File(descriptorooe_hdf5_path, 'r')
    term_full_ooe = descriptor_hdf5['term_roi_ooe'][:]
    term_full_counts = descriptor_hdf5['term_roi_counts'][:]
    descriptor_hdf5.close()

    #total number of terms in all regions
    total_term_counts = np.sum(term_full_counts)
    #split on cluster and find significant terms
    print('Determining significance for terms in each region...', flush=True)
    func = partial(find_significance_for_region, total_term_counts)
    with closing(Pool(num_cpus)) as pool:
        sig_pa = pool.map(func, np.arange(len(term_full_counts)))
    pool.join()
    sig_pa = np.array(sig_pa)

    #keep only ooe values that are also significant
    term_sig_ooe = np.where(sig_pa, np.log2(term_full_ooe), 0)
    term_sig_counts = np.where(sig_pa, term_full_counts, 0)

    #save the terms and the term counts for each region
    out_f = h5py.File(out_path, 'w')
    dset = out_f.create_dataset('sig_term_roi_counts', data=term_sig_counts)
    dset = out_f.create_dataset('sig_term_roi_log2ooe', data=term_sig_ooe)
    out_f.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
