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


dhs_df = pd.read_pickle('../data/data_240.pkl')

dhs_arr = dhs_df.to_numpy()

dhs_dist = []

def combined_func(dhs_arr, index):
    distr = []
    print(f"start process {index}")
    for row in range(index, len(dhs_arr), 100):
        if row==index:
            continue
        else:
            distr.append(jaccard(dhs_arr[index], dhs_arr[row]))
    print(f"end process {index}")
    return np.array(distr)
    


func = partial(combined_func, dhs_arr)
with closing(Pool(10)) as pool:
    results = pool.map(func, list(range(0, len(dhs_arr), 100)))
pool.join()

distr = np.concatenate(results, axis=0)

results = sorted(results, key = len)               
                   
pkl_file = open('../output/dhs-jaccard-240s-scores-every-100.pkl','wb')
pickle.dump(results, pkl_file)

fn = '../output/dhs-jaccard-240s-distr-every-100.npy'
np.save(fn, distr)

distr_plot = sns.displot(distr, height=5, aspect=3)
distr_plot.set(title='DHS Jaccard Distance (240S) Distribution')
plt.savefig('../output/plots/ST_' + today.strftime("%Y-%m-%d") + '_DHS-Jaccard-240-sample-distribution-every-100.png')



# with open('../output/DHS_distribution.txt', 'w'):
    