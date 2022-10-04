import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from functools import partial
from multiprocessing import Pool
from contextlib import closing
import sys
import random
from scipy.stats import percentileofscore
import os

#find the cosine simliarity for 1 region compared to all regions
def cosine_sim_per_ind(ooe_arr, sample_inds, tmp_path, ind):
    cosine_sim_arr = cosine_similarity(ooe_arr[[ind], :], ooe_arr)[0]
    if ind in sample_inds:
        np.save(os.path.join(tmp_path, 'sample_'+str(ind)+'_cosine_arr.npy'), cosine_sim_arr)
    cosine_sim_arr[ind] = -1
    cosine_sim_arr = (-np.sort(-cosine_sim_arr))[:100]
    return cosine_sim_arr

#find the cosine simliarity for 1 region compared to all regions in terms of region inds
def cosine_sim_regions_per_ind(ooe_arr, sample_inds, tmp_path, ind):
    cosine_sim_arr = cosine_similarity(ooe_arr[[ind], :], ooe_arr)[0]
    cosine_sim_arr[ind] = -1
    cosine_sim_inds_arr = (np.argsort(-cosine_sim_arr))[:100]
    return cosine_sim_inds_arr

#find the percentiles of top 100 cosine values for each region
def percentile_calc(cosine_samples, sample_top_cosine):
    #sample_top_cosine = top_100_cosine_vals[ind, :]
    sample_percentiles = []
    for score in sample_top_cosine:
        sample_percentiles.append(percentileofscore(cosine_samples, score))
    return sample_percentiles

def main(ooe_path, num_cpus, cosine_path, tmp_path):
    #load the data
    print('Loading data...', flush=True)
    ooe_arr = np.load(ooe_path)
    
    sample_inds = np.array(random.sample(range(len(ooe_arr)), k=100))

    #split on cluster and calculate cosine similarity for all regions
    print('Calculating cosine similarity...', flush=True)
    func1 = partial(cosine_sim_per_ind, ooe_arr, sample_inds, tmp_path)
    func2 = partial(cosine_sim_regions_per_ind, ooe_arr, sample_inds, tmp_path)
    with closing(Pool(num_cpus)) as pool:
        top_100_cosine_vals = pool.map(func1, np.arange(len(ooe_arr)))
        top_100_similar_inds = pool.map(func2, np.arange(len(ooe_arr)))
    pool.join()
    top_100_cosine_vals = np.array(top_100_cosine_vals)
    top_100_similar_inds = np.array(top_100_similar_inds)

    #sample 100 regions in order to calculate percentiles
    print('Sampling...')
    #sample_inds = np.array(random.sample(range(len(ooe_arr)), k=100))
    cosine_samples = np.array([])
    for i in sample_inds:
        cosine_samples = np.append(cosine_samples, np.load(os.path.join(tmp_path, 'sample_'+str(i)+'_cosine_arr.npy')))
    #cosine_samples = np.ndarray.flatten(cosine_m[sample_inds, :])

    #find the top 100 most similar regions for each region 
    #print('Finding top 100 regions...', flush=True)
    #cosine_m[np.arange(len(ooe_arr)), np.arange(len(ooe_arr))] = -1
    
    #top_100_cosine_vals = -np.sort(-cosine_m, axis=1)[:, :100]
    #top_100_simlar_inds = np.argsort(-cosine_m, axis=1)[:, :100]
    #top_100_cosine_vals = []
    #top_100_similar_inds = []

    #while len(cosine_m) > 0:
     #   row = cosine_m[0, :]
      #  top_100_cosine_vals.append((-np.sort(-row))[:100])
       # top_100_similar_inds.append((np.argsort(-row))[:100])
       # cosine_m = cosine_m[1:, :]
    
    #top_100_cosine_vals = np.array(top_100_cosine_vals)
    #top_100_similar_inds = np.array(top_100_similar_inds)

    #del cosine_m

    #generate percentiles of cosine similarities for top 100 most similar regions
    print('Calculating percentiles...', flush=True)
    func3 = partial(percentile_calc, cosine_samples)
    with closing(Pool(num_cpus)) as pool:
        top_100_percentiles = pool.map(func3, top_100_cosine_vals[np.arange(len(ooe_arr)), :])
    pool.join()
    top_100_percentiles = np.array(top_100_percentiles)

    #save top 100 of cosine similarities, similar regions, and percentiles for each region
    np.savez(cosine_path, top_100_cosine_vals=top_100_cosine_vals, top_100_similar_regions=top_100_similar_inds, top_100_percentile_vals=top_100_percentiles)

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4])
