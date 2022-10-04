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


chromatin_state_arr = np.load('../data/epilogos_state_data_18_1000_arr.npy')

chromatin_state_dist = []

def combined_func(chromatin_state_arr, index):
    distr = []
    print(f"start process {index}")
    for row in range(index, len(chromatin_state_arr)):
        if row==index:
            continue
        else:
            distr.append(cosine(chromatin_state_arr[index], chromatin_state_arr[row]))
    print(f"end process {index}")
    return np.array(distr)
    


func = partial(combined_func, chromatin_state_arr)
with closing(Pool(10)) as pool:
    results = pool.map(func, list(range(0, len(chromatin_state_arr))))
pool.join()

distr = np.concatenate(results, axis=0)

results = sorted(results, key = len)               
                   
pkl_file = open('../output/chromatin_state_cosine_distr_1000bp_every1000.pkl','wb')
pickle.dump(results, pkl_file)

fn = '../output/chromatin_state_cosine_distr_200bp_every1000.npy'
np.save(fn, distr)

distr_plot = sns.displot(distr, height=5, aspect=3)
distr_plot.set(title='Every 100th 200bp bin chromatin state cosine distribution')
plt.savefig('../output/plots/ST_' + today.strftime("%Y-%m-%d") + '_chromatin-state-cosine-distribution-every-1000.png')



# with open('../output/DHS_distribution.txt', 'w'):
    