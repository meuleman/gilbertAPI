import numpy as np
from functools import partial
from multiprocessing import Pool
from contextlib import closing
import sys
import h5py
from collections import Counter
import re
from scipy.stats import hypergeom


def find_significance_for_region(region):
    #for a single region, check if each term that occurs is significant and return mask array
    region_counts = gene_set_full_counts[region]
    a = region_counts
    b = np.array([101]*len(region_counts))
    c = c_sum_all
    d = np.array([101*len(gene_set_full_counts)]*len(region_counts))
    conting = np.array([[a, b], 
                        [c, d]])
    pvals = fishers_test_greater(conting)
    #bonferroni correction
    sig_pa = pvals < (0.05 / len(region_counts))
    return sig_pa

def fishers_test_greater(table):
    #code modified from https://github.com/scipy/scipy/blob/v1.8.1/scipy/stats/_stats_py.py#L4120-L4372
    #scipy stats fisher's exact test, but so that we can perform multiple tests at once
    c = np.asarray(table, dtype=np.int64)
    
    if np.any(c < 0):
        raise ValueError("All values in `table` must be nonnegative.")
        
    n1 = c[0, 0, :] + c[0, 1, :]
    n2 = c[1, 0, :] + c[1, 1, :]
    n = c[0, 0, :] + c[1, 0, :]
    
    pvalue = hypergeom.cdf(c[0, 1, :], n1 + n2, n1, c[0, 1, :] + c[1, 1, :])
    pvalue = np.where(pvalue > 1, 1, pvalue)
    
    return pvalue

def main(genesetooe_hdf5_path, out_path, num_cpus):
    #load data
    global gene_set_full_counts
    geneset_hdf5 = h5py.File(genesetooe_hdf5_path, 'r')
    gene_set_full_ooe = geneset_hdf5['gene_set_roi_ooe'][:]
    gene_set_full_counts = geneset_hdf5['gene_set_roi_counts'][:]
    geneset_hdf5.close()

    global c_sum_all
    c_sum_all = np.sum(gene_set_full_counts, axis=0)

    #split on cluster and find significant terms
    print('Determining significance for terms in each region...', flush=True)
    func = partial(find_significance_for_region)
    with closing(Pool(num_cpus)) as pool:
        sig_pa = pool.map(func, np.arange(len(gene_set_full_counts)))
    pool.join()
    sig_pa = np.array(sig_pa)

    #keep only ooe values that are also significant
    gene_set_sig_ooe = np.where(sig_pa, np.log2(gene_set_full_ooe), 0)
    gene_set_sig_counts = np.where(sig_pa, gene_set_full_counts, 0)

    #save the terms and the term counts for each region
    out_f = h5py.File(out_path, 'w')
    dset = out_f.create_dataset('sig_gene_set_roi_counts', data=gene_set_sig_counts)
    dset = out_f.create_dataset('sig_gene_set_roi_log2ooe', data=gene_set_sig_ooe)
    out_f.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
