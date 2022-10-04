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

columns = ['seqname', 'start', 'end', '?']
bins_25kb = pd.read_csv('../data/Altius_100k_epilogos_hg19_WM20210518.bed', delimiter='\t', names=columns)

tf_df = pd.read_pickle('../data/TF-df-hg38.pkl')
print('read in transcription factor df')

# Calculate expected values
# Number of regions:
num_regions = len(bins_25kb)
# Total number of DHSs per component:
expected = tf_df['TF motif cluster/group'].value_counts()/num_regions

# enrichment = enrichment.T
# enrichment = enrichment.set_index(list(range(len(bins_25kb))))

enrichment = pd.DataFrame()

def combined_func(tf_df, bins_25kb, index):
    print(f"start process {index}")
    
    # Look for all TFs in 25kb bins
    region_start = bins_25kb.loc[index].start
    region_end = bins_25kb.loc[index].end
    chromosome = bins_25kb.loc[index].seqname
    tf_chr = dhs_metadata[dhs_metadata['chromosome'] == chromosome]
    tf_sub = tf_chr.loc[(tf_chr['start'] > region_start) & (tf_chr['end'] < region_end)]
    del tf_chr
    tf_motif_cluster_count = tf_sub['TF motif cluster/group'].value_counts()
    del tf_sub
    enrichment = tf_motif_cluster_count.divide(expected)
    print(f"end process {index}")
    return enrichment
    


func = partial(combined_func, tf_df, bins_25kb)
with closing(Pool(10)) as pool:
    results = pool.map(func, list(range(0, len(bins_25kb), 1000)))
pool.join()

enrichment = pd.concat(results, axis=1)

enrichment = np.log2(enrichment)             
                   
pkl_file = open('../output/tf_enrichment_25kb.pkl','wb')
pickle.dump(results, pkl_file)


# with open('../output/DHS_distribution.txt', 'w'):
    