import numpy as np
import sys
from scipy.stats import rankdata
from scipy.spatial.distance import pdist, squareform
import h5py


# calculates sigmoid value of rank
def sigmoid_func(ranked_arr, scaling_factor=5):

    # max rank
    max_rank = np.max(ranked_arr, axis=0)

    # scaling between [-scaling_factor, scaling_factor]
    x = (ranked_arr / (max_rank / (scaling_factor * 2))) - scaling_factor
    sig_arr = 1 / (1 + np.exp(-x))
    return sig_arr


# top 1000 regions and values
def find_top_regions(sim_m):
    top_regions_all = []
    top_vals_all = []
    for i in range(len(sim_m)):
        top_regions = np.argsort(-sim_m[i])[:1001]
        top_regions = top_regions[top_regions != i][:1000]
        top_regions_all.append(top_regions)
        top_vals_all.append(sim_m[i, top_regions])
    top_regions_all = np.array(top_regions_all)
    top_vals_all = np.array(top_vals_all)
    return top_regions_all, top_vals_all


def rank_regions(arr):
    # rank each region
    ranked_regions = rankdata(arr, axis=0)

    # find the number of regions with 0 values for each component
    comp_num_zero = len(arr) - np.count_nonzero(arr, axis=0)

    # adjusting/resetting ranks
    ranked_regions = ranked_regions - comp_num_zero
    ranked_regions = np.where(ranked_regions < 0, 0, ranked_regions)

    return ranked_regions


def scale_by_sigmoid(sigmoid_ranks, top_regions_all, top_vals_all):
    # max sigmoid value across all components for all regions
    sigmoid_ranks_max = np.max(sigmoid_ranks, axis=1)
    # region a max sigmoid value
    max_sigmoid_a = np.vstack([sigmoid_ranks_max]*1000).T
    # region b max sigmoid value
    max_sigmoid_b = sigmoid_ranks_max[top_regions_all]
    # average max for a and b
    max_sigmoid_to_scale = (max_sigmoid_a + max_sigmoid_b) / 2
    # scale similarity scores
    top_sigmoid_scaled_vals_all = top_vals_all * max_sigmoid_to_scale
    return top_sigmoid_scaled_vals_all


def main(val_path, num_cpus, out_path):
    # load the data
    print('Loading data...', flush=True)
    val_hdf5 = h5py.File(val_path, 'r')
    val_arr = val_hdf5['base_arr'][:]
    val_hdf5.close()

    # for epilogos scores, setting negative values to 0
    val_arr = np.where(val_arr < 0, 0, val_arr)

    print('Ranking...', flush=True)
    # rank regions
    ranked_regions = rank_regions(val_arr)

    # max rankings
    max_ranks = np.max(ranked_regions, axis=0)

    print('Sigmoid...', flush=True)
    # sigmoid of rank for all regions
    sigmoid_ranks = sigmoid_func(ranked_regions)

    print('Euclidean...', flush=True)
    # euclidean distance matrix
    distance_m = squareform(pdist(sigmoid_ranks, metric='euclidean'))

    # find the max distance
    max_distance = np.max(distance_m)

    print('Similarity...', flush=True)
    # sim matrix
    # scaling between 0 and 1 by the max observed value
    sim_m = 1 - (distance_m / max_distance)
    del distance_m

    # regions with 0 datapoints
    total_val_per_region = np.sum(val_arr, axis=1)
    regions_0_val = np.where(total_val_per_region == 0)[0]

    # setting similarity values in regions with no datapoints to 0
    sim_m[regions_0_val] = 0
    sim_m[:, regions_0_val] = 0

    print('Top regions...', flush=True)
    # find top 1000 regions and values
    top_regions_all, top_vals_all = find_top_regions(sim_m)

    print('Sigmoid scaling...', flush=True)
    # scaling similarity scores by average of max
    # sigmoid values in both regions
    top_sigmoid_scaled_vals_all = scale_by_sigmoid(sigmoid_ranks,
                                                   top_regions_all,
                                                   top_vals_all)

    print('Saving...', flush=True)
    # save top 100 of sigmoid euclidean similarities, similar regions, and
    # percentiles for each region
    out_f = h5py.File(out_path, 'w')
    out_f.create_dataset('top_1000_scaled_vals', data=top_vals_all)
    out_f.create_dataset('top_1000_similar_regions', data=top_regions_all)
    out_f.create_dataset('top_1000_sigmoid_scaled_vals',
                         data=top_sigmoid_scaled_vals_all)
    out_f.create_dataset('sigmoid_ranks', data=sigmoid_ranks)
    out_f.create_dataset('ranked_regions', data=ranked_regions)
    out_f.create_dataset('max_ranks', data=max_ranks)
    out_f.close()


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
