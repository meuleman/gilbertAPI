import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from functools import partial
from multiprocessing import Pool
from contextlib import closing
import sys
import random
from scipy.stats import rankdata
import os

#find the cosine simliarity for 1 region compared to all regions in terms of region inds
def cosine_sim_regions_per_ind(ooe_arr, ind):
    cosine_sim_arr = cosine_similarity(ooe_arr[[ind], :], ooe_arr)[0]
    cosine_sim_arr[ind] = -1
    cosine_sim_inds_arr = (np.argsort(-cosine_sim_arr))[:1000]
    return cosine_sim_inds_arr

#find the cosine simliarity for 1 region compared to all identified regions
def cosine_sim_per_ind(ooe_arr, top_1000_similar_inds, ind):
    cosine_sim_arr = cosine_similarity(ooe_arr[[ind], :], ooe_arr[top_1000_similar_inds[ind], :])[0]
    return cosine_sim_arr

def main(ooe_path, num_cpus, cosine_path):
    #load the data
    print('Loading data...', flush=True)
    ooe_arr = np.load(ooe_path)
    
    #split on cluster and find similar regions based on cosine similarity
    print('Finding similar regions...', flush=True)
    func = partial(cosine_sim_regions_per_ind, ooe_arr)
    with closing(Pool(num_cpus)) as pool:
        top_1000_similar_inds = pool.map(func, np.arange(len(ooe_arr)))
    pool.join()
    top_1000_similar_inds = np.array(top_1000_similar_inds)

    #split on cluster and find cosine similarity values for identified regions
    print('Finding cosine similarity values...', flush=True)
    func = partial(cosine_sim_per_ind, ooe_arr, top_1000_similar_inds)
    with closing(Pool(num_cpus)) as pool:
        top_1000_cosine_vals = pool.map(func, np.arange(len(ooe_arr)))
    pool.join()
    top_1000_cosine_vals = np.array(top_1000_cosine_vals)

    #ranking cosine values and setting rank score of cosine values=0 to 0
    print('Ranking Cosine Scores...', flush=True)
    top_1000_rank_vals = (1000 - (rankdata(-top_1000_cosine_vals, method='max', axis=1) - 1)) / 1000
    top_1000_rank_vals = np.where(top_1000_cosine_vals!=0, top_1000_rank_vals, 0)

    #save top 100 of cosine similarities, similar regions, and percentiles for each region
    np.savez(cosine_path, top_1000_cosine_vals=top_1000_cosine_vals, top_1000_similar_regions=top_1000_similar_inds, top_1000_rank_vals=top_1000_rank_vals)

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
