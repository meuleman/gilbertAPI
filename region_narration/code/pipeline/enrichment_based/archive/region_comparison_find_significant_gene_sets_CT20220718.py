import numpy as np
from functools import partial
from multiprocessing import Pool
from contextlib import closing
import sys
import h5py
from collections import Counter
import re
from scipy.stats import fisher_exact

def find_significance_for_region(total_gene_set_counts, region):
    #for a single region, check if each term that occurs is significant and return mask array
    region_counts = gene_set_full_counts[region]
    gene_set_inds_to_test = np.where(region_counts!=0)[0]
    sig_pa = np.array([False]*len(region_counts))
    for gene_set in gene_set_inds_to_test:
        gene_set_count_in_region = region_counts[gene_set]
        not_gene_set_count_in_region = np.sum(region_counts) - gene_set_count_in_region
        gene_set_count_in_not_region = np.sum(gene_set_full_counts[:, gene_set]) - gene_set_count_in_region
        not_gene_set_count_in_not_region = total_gene_set_counts - gene_set_count_in_region - not_gene_set_count_in_region - gene_set_count_in_not_region
        #calculating p value from fisher exact test
        conting = [[gene_set_count_in_region, not_gene_set_count_in_region], 
                   [gene_set_count_in_not_region, not_gene_set_count_in_not_region]]
        odds, pval = fisher_exact(conting, alternative='greater')
        #bonferroni correction
        if pval < (0.05 / len(gene_set_inds_to_test)):
            sig_pa[gene_set] = True
    return sig_pa

def main(genesetooe_hdf5_path, out_path, num_cpus):
    #load data
    global gene_set_full_counts
    geneset_hdf5 = h5py.File(genesetooe_hdf5_path, 'r')
    gene_set_full_ooe = geneset_hdf5['gene_set_roi_ooe'][:]
    gene_set_full_counts = geneset_hdf5['gene_set_roi_counts'][:]
    geneset_hdf5.close()

    #total number of terms in all regions
    total_gene_set_counts = np.sum(gene_set_full_counts)
    #split on cluster and find significant terms
    print('Determining significance for terms in each region...', flush=True)
    func = partial(find_significance_for_region, total_gene_set_counts)
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
