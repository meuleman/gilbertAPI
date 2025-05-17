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


def main():
    datadir = '/home/ctrader/region_similarity/data/rank_based/'
    in_f = datadir + 'regions_25kb.search_space_top_1000_no_overlap.mapped_extended_domain_gene_h_c2_genesets.bin.CT20230531.npy'
    search_space_geneset_bin = np.load(in_f)

    clustering_method='jaccard'
    print('Start', flush=True)
    search_space_geneset_bin_jac = pdist(search_space_geneset_bin, metric=clustering_method)
    print('Jaccard', flush=True)
    search_space_geneset_z_cluster = hierarchy.average(search_space_geneset_bin_jac)
    print('average', flush=True)
    search_space_geneset_z_cluster_opt = hierarchy.optimal_leaf_ordering(search_space_geneset_z_cluster, search_space_geneset_bin_jac)
    print('optimal ordering', flush=True)

    search_space_geneset_z_cluster_opt_leaves = hierarchy.leaves_list(search_space_geneset_z_cluster_opt)
    print('leaves', flush=True)

    # saving
    print('Saving', flush=True)
    datadir = '/home/ctrader/region_similarity/data/rank_based/'
    out_f = datadir + 'regions_25kb.search_space_top_1000_no_overlap.mapped_extended_domain_gene_h_c2_genesets.bin.regions_jac_ave_clusters.CT20230531.npz'

    np.savez(out_f, 
             search_space_geneset_bin_jac=search_space_geneset_bin_jac, 
             search_space_geneset_z_cluster=search_space_geneset_z_cluster,
             search_space_geneset_z_cluster_opt=search_space_geneset_z_cluster_opt,
             search_space_geneset_z_cluster_opt_leaves=search_space_geneset_z_cluster_opt_leaves)


if __name__ == "__main__":
    main()
