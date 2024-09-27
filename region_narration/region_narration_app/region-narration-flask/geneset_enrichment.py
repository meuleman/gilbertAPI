from flask import Blueprint, request
import numpy as np
import json
import h5py
import os
import time
from scipy.stats import hypergeom


# Create a Blueprint
geneset_enrichment_bp = Blueprint('geneset_enrichment_bp', __name__)

# Path to binary matrix file
datadir="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/geneset_enrichment/20240918"
bin_path = os.path.join(datadir, "c5.go.bp.v2024.1.bin.h5")


"""
Loads the binary matrix data from h5 file
"""
def load_geneset_data(bin_path):
    # load binary matrix, gene names, and geneset names
    with h5py.File(bin_path, 'r') as f:
        gsXgene = f['bin'][:]
        genes = np.char.decode(f['genes'][:])
        genesets = np.char.decode(f['genesets'][:])

    return gsXgene, genes, genesets

# preload the genesetXgene binary matrix data
gsXgene, genes, genesets = load_geneset_data(bin_path)


"""
Performs hypergeometric test with provided list of genes to determine enrichment for each geneset.
Returns a list of p-values for each geneset.
"""
def perform_geneset_hypergeometric_test(genes_for_test, gsXgene=gsXgene, genes=genes, genesets=genesets):
    start = time.time()

    # take the unique genes (in case neighboring regions point to the same gene)
    genes_for_test = np.unique(genes_for_test)

    # perform hypergeometric test for each geneset
    M = len(genes)  # total number of genes
    N = len(genes_for_test)  # number of provided genes to test geneset enrichment with
    n = np.count_nonzero(gsXgene, axis=1)  # number of genes in each geneset
    k = np.count_nonzero(gsXgene[:, np.isin(genes, genes_for_test)], axis=1)  # number of the provided genes in each geneset

    # Perform the hypergeometric test
    p = hypergeom.sf(k-1, M, n, N)

    # # bonferroni correction
    # p = (p * len(genesets)).clip(max=1)

    print(f"Time to calculate enrichments for each geneset: {time.time() - start}", flush=True)
    return p


"""
Sorts the p-values and genesets and returns a list of genesets with p-values less than threshold
"""
def generate_return_data(p, genesets=genesets, threshold=0.05):
    start = time.time()

    # sort by p-value
    sort_order = np.argsort(p)
    p_sorted = p[sort_order]
    genesets_sorted = genesets[sort_order]

    return_list = [{"geneset": gs, "p": float(p)} for gs, p in zip(genesets_sorted, p_sorted) if p < threshold]

    print(f"Time to generate return data: {time.time() - start}", flush=True)
    return json.dumps(return_list)


"""
Collects the genes in a geneset
"""
def get_geneset_genes(geneset):
    geneset_genes = genes[gsXgene[genesets == geneset][0]]
    return geneset_genes.tolist()



"""
Endpoint for calculating enrichment for genesets given a list of genes
"""
@geneset_enrichment_bp.route("geneset_enrichment", methods=['POST'])
def geneset_enrichmenr():
    request_info = request.get_json()
    genes_for_test = request_info['genes']  # ["gene1", "gene2", ...]
    threshold = request_info.get('threshold', 0.05)  # default threshold is 0.05

    # test enrichment of genes 
    p = perform_geneset_hypergeometric_test(genes_for_test)

    # generate return data
    return_list = generate_return_data(p, threshold=threshold)

    return return_list


"""
Endpoint for inverse geneset search, where a geneset is provided and a list of genes is returned
"""
@geneset_enrichment_bp.route("geneset_to_genes", methods=['POST'])
def geneset_to_genes():
    request_info = request.get_json()
    geneset = request_info['geneset']

    # get genes in geneset
    geneset_genes = get_geneset_genes(geneset)

    return json.dumps(geneset_genes)