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

today = date.today()


def main():
    datadir = '/home/ctrader/region_similarity/data/rank_based/'
    in_f = datadir + 'h_c2_c5GO_genesets.GREATv4_genes.bin.filtered.CT20230601.npz'
    gsXg_bin = np.load(in_f)['gsXg_bin']

    clustering_method='jaccard'
    print('Start', flush=True)

    gsXg_bin_jac = pdist(gsXg_bin, metric=clustering_method)
    print('Jaccard', flush=True)

    gsXg_bin_jac_z_cluster = hierarchy.average(gsXg_bin_jac)
    print('average', flush=True)

    gsXg_bin_jac_z_cluster_opt = hierarchy.optimal_leaf_ordering(gsXg_bin_jac_z_cluster, gsXg_bin_jac)
    print('optimal ordering', flush=True)

    gsXg_bin_jac_z_cluster_opt_leaves = hierarchy.leaves_list(gsXg_bin_jac_z_cluster_opt)
    print('leaves', flush=True)

    # saving
    print('Saving', flush=True)
    out_f = datadir + 'h_c2_c5GO_genesets.GREATv4_genes.bin.filtered.jac_ave_clusters.CT20230601.npz'
    np.savez(out_f,
             gsXg_bin_jac=gsXg_bin_jac, 
             gsXg_bin_jac_z_cluster=gsXg_bin_jac_z_cluster,
             gsXg_bin_jac_z_cluster_opt=gsXg_bin_jac_z_cluster_opt,
             gsXg_bin_jac_z_cluster_opt_leaves=gsXg_bin_jac_z_cluster_opt_leaves)


if __name__ == "__main__":
    main()
