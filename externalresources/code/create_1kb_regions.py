import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
from datetime import date
import re
from tqdm import tqdm
import time
from pybedtools import BedTool
import random
import math
import pickle
import seaborn as sns
from scipy.stats import percentileofscore
from scipy.spatial.distance import jaccard
from multiprocessing import Pool
from contextlib import closing
from functools import partial

today = date.today()

bin_cols = ['seqname', 'start', 'end', '?']
bins_25kb = pd.read_csv('../data/Altius_100k_epilogos_hg19_WM20210518.bed', delimiter='\t', names=bin_cols)
print(bins_25kb)

bins_1kb = pd.DataFrame()
for i in bins_25kb.index:
    temp_bins = pd.DataFrame()
    seqname = bins_25kb.loc[i]['seqname']
    start = bins_25kb.loc[i]['start']
    for j in range(25):
        index = i*25 +j
        temp_bins[index] = ([seqname, start, start+1000])
        start = start + 1000
    bins_1kb = pd.concat([bins_1kb, temp_bins.T], axis=0)
    if i % 100 == 0:
        print(i)
print(bins_1kb)
bins_1kb.to_pickle('../data/bins_1kb.pkl')
    