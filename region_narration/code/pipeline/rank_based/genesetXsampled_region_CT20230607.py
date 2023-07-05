import numpy as np
import pandas as pd
import scipy.sparse as sparse
import itertools
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from datetime import date
import re
# from tqdm import tqdm
import time
# from pybedtools import BedTool
import random
import math
import pickle
import seaborn as sns
from scipy.stats import percentileofscore, pearsonr, rankdata, describe, mode, fisher_exact
import scipy as sp
import scipy.stats as stats
from scipy.spatial.distance import jaccard, euclidean, cosine, pdist
# import functools
# import operator
import gzip
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances
from sklearn.metrics import pairwise_distances
import sys
import os
from functools import partial
from contextlib import closing
from multiprocessing import Pool
from sklearn.decomposition import NMF, non_negative_factorization
from collections import Counter
from scipy.stats import hypergeom
from scipy.stats.mstats import gmean
from scipy.cluster import hierarchy
from subprocess import Popen, PIPE
import requests
import json

today = date.today()

def production_post_request(coords, scale, num_sim_regions):
    prod_url = 'https://explore.altius.org:5001/narration'
    post_body = {
        'location': coords,
        'factors': '[' + ','.join(np.arange(34).astype(str)) + ']',
        'roiURL': 'False',
        'scale': scale,
        'queryFactorThresh': '0.5',
        'numSimilarRegions': str(num_sim_regions),
        'regionMethod': 'sliding_window'
    }
    response = requests.post(prod_url, json=post_body)
    return response

def query_sim_hypergeom(region_inds, N, M, n, geneset_order, search_space_geneset_bin):
    k = np.sum(search_space_geneset_bin[region_inds], axis=0)
    inds_to_test = np.where(k)[0]

    test = hypergeom.sf(k[inds_to_test], M, n[inds_to_test], N)
    tested_genesets = geneset_order[inds_to_test]
    test = test * len(geneset_order)  # bonferroni correction
    test = np.where(test > 1, 1, test)
    return tested_genesets, inds_to_test, test


def main(num_regions_to_sample, scale, num_sim_regions):
    datadir = '/home/ctrader/region_similarity/data/rank_based/'
    in_f = datadir + 'regions_{}kb.search_space_top_1000_no_overlap.mapped_extended_domain_GREATv4_gene_c5GOBP_genesets.bin.CT20230607.npy'.format(scale)
    search_space_geneset_bin = np.load(in_f)

    datadir = '/home/ctrader/region_similarity/data/rank_based/'
    c5gobp_f = datadir + 'MSigDB/c5.go.bp.v2023.1.Hs.symbols.gmt'

    # build dictionary of genesets mapping to a list of genes
    geneset2genes = {}
    for gs_f in [c5gobp_f]:  # hallmark_f, c2_f, c5go_f
        f = open(gs_f)
        for line in f.readlines():
            line = line.strip().split('\t')
            gs = line[0]
            genes = line[2:]
            geneset2genes[gs] = genes
        f.close()
    geneset_order = np.unique(list(geneset2genes.keys()))

    # query region coordinates
    datadir = '/home/ctrader/region_similarity/data/rank_based/'
    ss_f = datadir + 'regions_{}kb_CT20230427.DHS.factornormNMF.sum_enrich_corank.epilogos.statesum_rank.percentiles.limited_search_space_top_1000_no_overlap.npz'.format(scale)
    ss_coords = np.load(ss_f)['coordinates']

    # for hypergeometeric test
    gs_region_counts = np.sum(search_space_geneset_bin, axis=0)
    num_regions = search_space_geneset_bin.shape[0]

    # geneset x number of sampled regions
    geneset_sampled_region_nlog10p = np.zeros((len(geneset_order), num_regions_to_sample), dtype=float)
    # regions to sample
    sampled_region_inds = np.random.randint(0, num_regions, num_regions_to_sample)
    for sample_ind, sample_i in zip(sampled_region_inds, range(num_regions_to_sample)):
        c = ss_coords[sample_ind]
        query_coords = c[0] + ':' + c[1] + '-' + c[2]

        # get similar regions
        response = production_post_request(query_coords, scale, num_sim_regions)
        query_full_data = json.loads(response.text)

        # find best detail level
        prod_ranks = []
        for i in range(len(query_full_data)):
            prod_ranks.append(query_full_data[i][0]['prod_rank'])
        dl = np.argmax(prod_ranks) + 1

        query_sim_region_inds = []
        query_sim_region_inds.append(sample_ind)
        for i in range(1, len(query_full_data[0])):
            query_sim_region_inds.append(query_full_data[dl - 1][i]['region_ind'])

        # test for geneset significance
        query_tested_gs, query_tested_gs_inds, query_tested_gs_pval = query_sim_hypergeom(
            query_sim_region_inds,
            N=(num_sim_regions + 1),
            M=num_regions,
            n=gs_region_counts,
            geneset_order=geneset_order,
            search_space_geneset_bin=search_space_geneset_bin
        )
        # significant genesets for query region
        sig_inds = np.where(query_tested_gs_pval < 0.05)[0]
        query_sig_gs_inds = query_tested_gs_inds[sig_inds]
        query_sig_gs_pval = query_tested_gs_pval[sig_inds]
        # save -log10 pval for each significant geneset
        geneset_sampled_region_nlog10p[query_sig_gs_inds, sample_i] = -np.log10(query_sig_gs_pval)

    datadir = '/home/ctrader/region_similarity/data/rank_based/'
    out_f = datadir + 'regions_{}kb_search_space_c5GOBP_geneset_sampling.{}_regions.{}_sim_regions.CT20230607.npz'.format(scale, num_regions_to_sample, num_sim_regions)
    np.savez(out_f, geneset_sampled_region_nlog10p=geneset_sampled_region_nlog10p, geneset_order=geneset_order, sampled_region_inds=sampled_region_inds)

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
