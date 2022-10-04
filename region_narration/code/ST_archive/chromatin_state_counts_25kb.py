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
from scipy.spatial.distance import jaccard, cosine
from multiprocessing import Pool
from contextlib import closing
from functools import partial

today = date.today()

columns = ['seqname', 'start', 'end', '?']
bins_25kb = pd.read_csv('../data/Altius_100k_epilogos_hg19_WM20210518.bed', delimiter='\t', names=columns)

chromatin_state_data = pd.read_pickle('../data/epilogos_state_data_18.pkl')

print('files read')

print(chromatin_state_data.shape)

def combined_func(chromatin_state_data, index):
    print(f"start process {index}")
    region_start = bins_25kb.loc[row].start
    region_end = bins_25kb.loc[row].end
    chromosome = bins_25kb.loc[row].seqname
    chromatin_state_sub = chromatin_state_data[chromatin_state_data['seqname'] == chromosome]
    chromatin_state_sub = chromatin_state_sub.loc[(chromatin_state_sub['start'] > region_start) & (chromatin_state_sub['end'] < region_end)]
    print(f"end process {index}")
    return chromatin_state_sub.iloc[:, 3:].T.apply(pd.Series.value_counts, axis=0).sum(axis=1) 


func = partial(combined_func, chromatin_state_data)
with closing(Pool(8)) as pool:
    results = pool.map(func, list(range(0, len(chromatin_state_data), 100000)))
pool.join()

num_per_state_per_25kb = pd.concat(results, axis=1)

                   
pkl_file = open('../output/chromatin_state_count_every100000.pkl','wb')
pickle.dump(results, pkl_file)

# with open('../output/DHS_distribution.txt', 'w'):
    