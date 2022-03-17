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

dhs_metadata = pd.read_pickle('../data/dhs_metadata.pkl')

# Calculate expected values
# Number of regions:
num_regions = len(bins_25kb)
# Total number of DHSs per component:
dhs_per_comp = dhs_metadata.component.value_counts()
# Total number of DHS per component divided by Number of Regions
expected = dhs_per_comp/num_regions

# enrichment = enrichment.T
# enrichment = enrichment.set_index(list(range(len(bins_25kb))))

enrichment = pd.DataFrame()

def combined_func(dhs_metadata, bins_25kb, index):
    print(f"start process {index}")
    # Look for all DHSs in 25kb bins
    region_start = bins_25kb.loc[index].start
    region_end = bins_25kb.loc[index].end
    chromosome = bins_25kb.loc[index].seqname
    dhs_metadata_sub = dhs_metadata[dhs_metadata['seqname'] == chromosome]
    dhs_metadata_sub = dhs_metadata_sub.loc[(dhs_metadata_sub['start'] > region_start) & (dhs_metadata_sub['end'] < region_end)]
    component_count = dhs_metadata_sub['component'].value_counts()
    enrichment = component_count.divide(expected, fill_value=0)
    print(f"end process {index}")
    return enrichment
    


func = partial(combined_func, dhs_metadata, bins_25kb)
with closing(Pool(10)) as pool:
    results = pool.map(func, list(range(0, len(bins_25kb))))
pool.join()

enrichment = pd.concat(results, axis=1)

enrichment = np.log2(enrichment)             
                   
pkl_file = open('../output/dhs_enrichment_25kb.pkl','wb')
pickle.dump(results, pkl_file)


# with open('../output/DHS_distribution.txt', 'w'):
    