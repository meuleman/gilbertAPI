import h5py
import numpy as np
import time
import os
from scipy.stats import hypergeom
from subprocess import Popen, PIPE
import json


def fetch_indices_scores(
    order, 
    chromosome, 
    factor_index, 
    dataset_name, 
    file_path
):
    """
    Fetches the order 14 positional indices for a given order, chromosome, factor index, 
    and dataset name. Returns an empty array if the indices are not found.
    """
    filter_key = f"{order}.{chromosome}.{dataset_name}.{factor_index}"
    with h5py.File(file_path, 'r') as f:
        try:
            indices_scores = f[filter_key][:]
        except:
            indices_scores = np.array([], dtype="int32")
    return indices_scores


def preload_preview_indices(index_path, max_order=11):
    """
    Function to preload all indices to be used for filter previews up through the max order
    """
    preload_time = time.time()
    preloaded_indices = {}
    
    with h5py.File(index_path, 'r') as f:

        # find all the datasets less than or equal to the max order
        datasets = list(f.keys())
        datasets = [d for d in datasets if int(d.split('.')[0]) <= max_order]
        
        # preload indices for each dataset
        for d in datasets:

            # create chromosome independent key
            order, chrm, dataset_name, factor_index = d.split('.')
            key = f"{order}.{dataset_name}.{factor_index}"

            # all factor instances across chromosomes are saved together in a list
            if key not in preloaded_indices:
                preloaded_indices[key] = []
            
            # add to factorXorder list
            preloaded_indices[key].append({
                "chromosome": chrm,
                "order": int(order), 
                "index": int(factor_index), 
                "dataset_name": dataset_name, 
                "indices": f[d][:, 0].astype("int32"), 
                "scores": f[d][:, 1].astype("float32")
            })

    # save the order through which indices have been preloaded
    preloaded_indices['max_order'] = max_order
            
    print(f"Time to preload indices for filter previews: {time.time() - preload_time}", flush=True)
    return preloaded_indices


def hilbert_pos_to_order(pos, frm, to):
    """
    Converts an array of hilbert positions from one order to another
    """
    return np.floor(pos / pow(4.0, frm - to)).astype(int)


def tabix_query(query, path):
    """
    Queries a tabix file for provided genomic coordinates (chr:start-end)
    Return all the overlapping hits
    """
    # Use subprocess with tabix and return hits
    process = Popen(['tabix', '-f', path, query], stdout=PIPE)
    
    hits = []
    for line in process.stdout:
        hits.append(str(np.char.decode(line)).strip().split('\t'))
    return hits


def load_geneset_data(bin_path):
    """
    Loads the binary matrix data from h5 file
    """
    # load binary matrix, gene names, and geneset names
    with h5py.File(bin_path, 'r') as f:
        gsXgene = f['bin'][:]
        genes = np.char.decode(f['genes'][:])
        genesets = np.char.decode(f['genesets'][:])

    return gsXgene, genes, genesets


def perform_geneset_hypergeometric_test(genes_for_test, gsXgene, genes, genesets, membership=False, bonferroni=False):
    """
    Performs hypergeometric test with provided list of genes to determine enrichment for each geneset.
    Returns a list of p-values for each geneset. If membership is True, returns a boolean array of 
    geneset membership for test genes.
    """
    start = time.time()

    # take the unique genes (in case neighboring regions point to the same gene)
    genes_for_test = np.unique(genes_for_test)

    # perform hypergeometric test for each geneset
    M = len(genes)  # total number of genes
    N = len(genes_for_test)  # number of provided genes to test geneset enrichment with
    n = np.count_nonzero(gsXgene, axis=1)  # number of genes in each geneset
    k = np.count_nonzero(gsXgene[:, np.isin(genes, genes_for_test)], axis=1)  # number of the provided genes in each geneset

    
    if membership:
        # find the genesets in which genes have membership
        p = k > 0
    else:
        # Perform the hypergeometric test
        p = hypergeom.sf(k-1, M, n, N)

        # bonferroni correction
        p = (p * len(genesets)).clip(max=1) if bonferroni else p

    print(f"Time to calculate enrichments for each geneset: {time.time() - start}", flush=True)
    return p


def generate_geneset_return_data(p, genesets, gsXgene, genes, threshold=0.05):
    """
    Sorts the p-values and genesets and returns a list of genesets with p-values less than threshold
    """
    start = time.time()

    if p.dtype == bool:
        # return genesets in which there is membership
        return_list = [{"geneset": gs} for gs, p in zip(genesets, p) if p]

    else:
        # sort by p-value
        sort_order = np.argsort(p)
        p_sorted = p[sort_order]
        genesets_sorted = genesets[sort_order]

        return_list = [
            {
                "geneset": gs, 
                "p": float(p), 
                "genes": genes[gsXgene[i]].tolist()
            } for gs, p, i in zip(genesets_sorted, p_sorted, sort_order) if p < threshold
        ]

    print(f"Time to generate geneset return data: {time.time() - start}", flush=True)
    return json.dumps(return_list)