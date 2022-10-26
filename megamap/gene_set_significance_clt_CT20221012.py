import numpy as np
from functools import partial
from multiprocessing import Pool
from contextlib import closing
import sys
import h5py
import pickle
import scipy.stats as stats
import pandas as pd

def geneset_corrs_info(gene_sets, megamap_genes):
    num_genes_per_gs_mm = np.array([])
    gs_corrs_mean_all_gs = []
    for gs in gene_sets:
        #find genes
        genes_in_gs = c1c2hallmarkgeneset2gene[gs]
        #make sure in mm
        genes_in_gs_mm = megamap_genes[np.isin(megamap_genes, genes_in_gs)]
        #count genes
        num_genes_per_gs_mm = np.append(num_genes_per_gs_mm, len(genes_in_gs_mm))
        #get correlation values for genes
        gs_corrs = corrs_arr[:, np.isin(megamap_genes, genes_in_gs_mm)]
        #average correlation values
        gs_corrs_mean = np.mean(gs_corrs, axis=1)
        gs_corrs_mean_all_gs.append(gs_corrs_mean)
    gs_corrs_mean_all_gs = np.array(gs_corrs_mean_all_gs).T
    return gs_corrs_mean_all_gs, num_genes_per_gs_mm

def pval_calc(gs_corrs_mean_all_gs, num_genes_per_gs_mm):
    #mean for all experiments
    exper_corrs_mean = np.mean(corrs_arr, axis=1)
    #mean replicated for all gene sets
    exper_corrs_mean_all_gs = np.array([exper_corrs_mean]*len(num_genes_per_gs_mm)).T
    #standard deviation for all experiments
    exper_corrs_std = np.std(corrs_arr, axis=1)
    #std replicated for all gene sets
    exper_corrs_std_all_gs = np.array([exper_corrs_std]*len(num_genes_per_gs_mm)).T

    #corrected experiment std
    exper_corrs_std_all_gs_corrected = exper_corrs_std_all_gs / np.sqrt(num_genes_per_gs_mm) #std of population / sqrt of number of genes

    #z score for all experiments and gene sets
    corrs_z = (gs_corrs_mean_all_gs - exper_corrs_mean_all_gs) / exper_corrs_std_all_gs_corrected

    #pvalue for all experiments and gene sets (one-sided greater, bonferroni)
    corrs_p = stats.norm.sf(corrs_z) * len(num_genes_per_gs_mm) * len(corrs_arr)
    corrs_p = np.where(corrs_p > 1, 1, corrs_p)
    return corrs_p

def load_gs_dict(file):
    pickle_in = open(file, 'rb')
    gene2gs = pickle.load(pickle_in)
    pickle_in.close()
    return gene2gs

def reverse_dict(gene2gs):
    gs2gene = {}
    for gene in gene2gs.keys():
        gene_sets = gene2gs[gene]
        for gene_set in gene_sets:
            if gene_set not in gs2gene:
                gs2gene[gene_set] = np.array([gene])
            else:
                gs2gene[gene_set] = np.append(gs2gene[gene_set], np.array([gene]))
    return gs2gene

def main(corrs_file, gene2hallmark_file, gene2c1_file, gene2c2_file, out_file):
    #c1 gene set 
    gene2c1geneset = load_gs_dict(gene2c1_file)
    #c2 gene set 
    gene2c2geneset = load_gs_dict(gene2c2_file)
    #hallmark gene set
    gene2hallmarkgeneset = load_gs_dict(gene2hallmark_file)

    #reversing dictionaries
    hallmarkgeneset2gene = reverse_dict(gene2hallmarkgeneset)
    c1geneset2gene = reverse_dict(gene2c1geneset)
    c2geneset2gene = reverse_dict(gene2c2geneset)

    #combining dictionaries
    global c1c2hallmarkgeneset2gene
    c1c2hallmarkgeneset2gene = hallmarkgeneset2gene
    c1c2hallmarkgeneset2gene.update(c1geneset2gene)
    c1c2hallmarkgeneset2gene.update(c2geneset2gene)
    num_gs = len(c1c2hallmarkgeneset2gene)
    gene_sets = np.sort(list(c1c2hallmarkgeneset2gene.keys()))

    #experiment x gene correlations
    corrs = pd.read_feather(corrs_file)
    megamap_genes = corrs.columns[4:]
    global corrs_arr
    corrs_arr = np.abs(corrs[megamap_genes].to_numpy()) #absolute value

    #mean correlation score and number of genes for each gene set 
    gs_corrs_mean_all_gs, num_genes_per_gs_mm = geneset_corrs_info(gene_sets, megamap_genes)

    #p-value calculation
    corrs_p = pval_calc(gs_corrs_mean_all_gs, num_genes_per_gs_mm)
    
    #save
    out = h5py.File(out_file, 'w')
    dset = out.create_dataset('experiment_gs_pvalues', data=corrs_p)
    dset = out.create_dataset('gene_sets', data=np.char.encode(gene_sets, 'utf-8'))
    out.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

