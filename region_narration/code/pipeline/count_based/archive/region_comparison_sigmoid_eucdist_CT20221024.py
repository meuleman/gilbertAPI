import numpy as np
from functools import partial
from multiprocessing import Pool
from contextlib import closing
import sys
import random
from scipy.stats import rankdata
from scipy.spatial.distance import pdist, squareform
import os
import h5py

#ranking and scaling data
def scale_data(data, max_poss):
    data_to_rank = data[data!=0]
    max_obs = np.max(data_to_rank)
    rank_data = rankdata(data_to_rank, method='average')
    rank_data = rank_data/(len(data_to_rank))
    rank_data = rank_data * max_obs/max_poss
    
    original_data_shape = data.shape
    data = data.astype(float).flatten()
    np.put(data, np.where(data!=0)[0], rank_data)
    data = data.reshape(original_data_shape)
    return data

#calculates sigmoid value of rank
def sigmoid_func(arr, scaling_factor=5):
    num_regions = len(arr)
    x = (arr / (num_regions / (scaling_factor * 2))) - scaling_factor
    sig_arr = 1 / (1 + np.exp(-x))
    return sig_arr

#top 1000 regions and values
def find_top_regions(sim_m):
    top_regions_all = []
    top_vals_all = []
    for i in range(len(sim_m)):
        top_regions = np.argsort(-sim_m[i])[:1001]
        top_regions = top_regions[top_regions!=i][:1000]
        top_regions_all.append(top_regions)
        top_vals_all.append(sim_m[i, top_regions])
    top_regions_all = np.array(top_regions_all)
    top_vals_all = np.array(top_vals_all)
    return top_regions_all, top_vals_all

def main(count_path, num_cpus, out_path):
    #load the data
    print('Loading data...', flush=True)
    count_hdf5 = h5py.File(count_path, 'r')
    count_arr = count_hdf5['count_arr'][:]
    count_hdf5.close()

    print('Ranking...', flush=True)
    #rank each region based on counts
    ranked_counts = rankdata(count_arr, axis=0)

    #limit component ranks to just ragions with at least 1 datapoint of that component
    ranked_counts = np.where(count_arr==0, 0, ranked_counts)

    print('Sigmoid...', flush=True)
    #sigmoid of rank for all regions
    sigmoid_ranks = sigmoid_func(ranked_counts)

    print('Euclidean...', flush=True)
    #euclidean distance matrix
    distance_m = squareform(pdist(sigmoid_ranks, metric='euclidean'))

    #find the max distance
    max_distance = np.max(distance_m)

    print('Similarity...', flush=True)
    #sim matrix
    sim_m = max_distance - distance_m
    del distance_m

    #regions with 0 datapoints
    total_counts_per_region = np.sum(count_arr, axis=1)
    regions_0_counts = np.where(total_counts_per_region==0)[0]

    #setting similarity values in regions with no datapoints to 0
    sim_m[regions_0_counts] = 0
    sim_m[:, regions_0_counts] = 0

    print('Top regions...', flush=True)
    ####could split on cluster if takes too long
    #find top 1000 regions and values
    top_regions_all, top_vals_all = find_top_regions(sim_m)

    print('Scaling...', flush=True)
    #scale values
    top_scaled_vals = scale_data(top_vals_all, max_distance)

    print('Saving...', flush=True)
    #save top 100 of sigmoid euclidean similarities, similar regions, and percentiles for each region
    out_f = h5py.File(out_path, 'w')
    dset = out_f.create_dataset('top_1000_scaled_vals', data=top_scaled_vals)
    dset = out_f.create_dataset('top_1000_similar_regions', data=top_regions_all)
    dset = out_f.create_dataset('top_1000_sigmoideuc_vals', data=top_vals_all)
    out_f.close()

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
