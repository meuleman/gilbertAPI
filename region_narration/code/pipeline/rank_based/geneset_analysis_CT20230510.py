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
# from scipy.spatial.distance import jaccard, euclidean, cosine
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
import re
import h5py
import tabix
import upsetplot
from wordcloud import WordCloud
from sklearn.decomposition import NMF, non_negative_factorization
from collections import Counter
from scipy.stats import hypergeom
from scipy.stats.mstats import gmean
from scipy.cluster import hierarchy
from subprocess import Popen, PIPE

today = date.today()


def load_factors(per_f):
    hdf5_f = h5py.File(per_f, 'r')
    factors = np.char.decode(hdf5_f['factors'][:])
    color = np.char.decode(hdf5_f['colors'][:])
    hdf5_f.close()
    return factors, color


def tabix_query(filename, query):
    # Use subprocess with tabix and return hits
    process = Popen(['tabix', '-f', filename, query], stdout=PIPE)

    hits = []
    for line in process.stdout:
        hits.append(np.char.decode(line.strip().split()))
    hits = np.array(hits)
    return hits

# determine roi from query coordinates
def load_roi(genomic_location, link_path):

    # determine chromosome, start, and end from query
    query_start_end = genomic_location.split(':')[1].split('-')
    query_start = int(query_start_end[0])
    query_end = int(query_start_end[1])

    # retrieve possible rois
    possible_coords = tabix_query(link_path, genomic_location)

    # find closest roi
    start_dist = np.abs(possible_coords[:, 1].astype(int) - query_start)
    end_dist = np.abs(possible_coords[:, 2].astype(int) - query_end)
    query_region = possible_coords[np.argmin(start_dist + end_dist)]

    roi = int(query_region[3])
    roi_percents = query_region[4:].astype('float32')
    roi_coords = query_region[:3]

    return roi, roi_percents, roi_coords


def main():
    datadir = '/home/ctrader/region_similarity/data/rank_based/'
    DHS_per_f = datadir + 'regions_1kb_CT20230427.DHS.factornormNMF.sum_enrich_corank_percentile.hdf5'
    epi_per_f = datadir + 'regions_1kb_CT20230427.epilogos.statesum_rank_percentile.hdf5'

    # load region percentiles
    DHS_factors, DHS_color = load_factors(DHS_per_f)
    epi_factors, epi_color = load_factors(epi_per_f)

    factor_order = np.append(DHS_factors, epi_factors)





    datadir = '/home/ctrader/region_similarity/data/'
    hallmark_f =  datadir + 'h.all.v7.5.1.symbols.gmt'
    c1_f = datadir + 'c1.all.v7.5.1.symbols.gmt'
    c2_f = datadir + 'c2.all.v7.5.1.symbols.gmt'


    geneset2genes = {}
    for gs_f in [hallmark_f, c1_f, c2_f]:
        f = open(gs_f)
        for line in f.readlines():
            line = line.strip().split('\t')
            geneset2genes[line[0]] = line[2:]


    datadir = '/home/ctrader/region_similarity/data/'
    genes_f = datadir + 'gencode.v41.basic.annotation.parsed.CT20221007.txt'
    f = open(genes_f)

    gene2coords = {}
    for line in f.readlines():
        line = line.strip().split('\t')
        if line[4] == 'gene':
            gene = line[3]
            if gene not in gene2coords:
                gene2coords[gene] = np.array([line[0] + ':' + line[6] + '-' + line[7]])  # promoter bp
            else:
                gene2coords[gene] = np.append(gene2coords[gene], np.array([line[0] + ':' + line[6] + '-' + line[7]]))



    datadir = '/home/ctrader/region_similarity/data/rank_based/'
    out_f = datadir + 'single_coord_mapping_genes_in_gs_promoter_25kb_region_ranks.gencode.v41.CT20230510.npz'
    gs_expected_roi = np.load(out_f)['indicies']
    gs_expected_ranks = np.load(out_f)['ranks']
    gs_expected_coords = np.load(out_f)['coordinates']


    datadir = '/home/ctrader/region_similarity/data/rank_based/indexed/'
    scale = '25'
    local_tabix_path = datadir + 'hg38.' + scale + 'kb.regions.CT20230427.sorted.rownum4.factorranks5t38.bed.gz'



    geneset_tracker = []
    enrichment_tracker = []

    for geneset in list(geneset2genes.keys()):
        print(geneset, flush=True)
        geneset_tracker.append(geneset)
        genes = geneset2genes[geneset]

        # get coordinates for genes in geneset
        gs_true_coords = []
        for g in genes:
            if g in gene2coords:
                coords = gene2coords[g]
                if ((len(coords) == 1) & (coords[0].split(':')[0][-1] != 'M')):
                    gs_true_coords.append(coords[0])

        # get ranks
        gs_rois = []
        gs_ranks = []
        gs_coords = []
        for coord in gs_true_coords:
            roi, roi_ranks, roi_coords = load_roi(coord, local_tabix_path)
            gs_rois.append(roi)
            gs_ranks.append(roi_ranks)
            gs_coords.append(roi_coords)

        # no region occurring twice
        gs_rois, inds = np.unique(gs_rois, return_index=True)
        gs_ranks = np.array(gs_ranks)[inds]
        gs_coords = np.array(gs_coords)[inds]

        # get obs/exp
        obs = np.mean(gs_ranks, axis=0)
        exp = np.mean(gs_expected_ranks, axis=0)

        ooe = obs / exp
        log2ooe = np.log2(ooe)
        enrichment_tracker.append(log2ooe)



    geneset_tracker = np.array(geneset_tracker)
    enrichment_tracker = np.array(enrichment_tracker)


    datadir = '/home/ctrader/region_similarity/data/rank_based/'
    outf = datadir + 'all_geneset_25kb_regional_enrichment_analysis_CT20230510.npz'
    np.savez(outf, enrichments=enrichment_tracker, genesets=geneset_tracker, factor_order=factor_order)

if __name__ == "__main__":
    main()
