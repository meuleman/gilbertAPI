import numpy as np
import sys
import h5py
from scipy.stats import rankdata


def sigmoid_func(val_arr, scaling_factor=5, method='factor'):
    # max rank
    if method == 'factor':
        max_val = np.max(val_arr, axis=0)
        # to prevent nan
        max_val = np.where(max_val == 0, 1, max_val)
    elif method == 'full':
        max_val = np.max(val_arr)

    # scaling between [-scaling_factor, scaling_factor]
    x = (val_arr / (max_val / (scaling_factor * 2))) - scaling_factor
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


def main(base_path, out_path):
    # load the data
    hdf5 = h5py.File(base_path, 'r')
    data = hdf5['base_arr'][:]
    hdf5.close()

    # Priors
    # max for each region
    data_max = np.max(data, axis=1)
    # rank each region based on max
    max_rank = rank_regions(data_max)
    # sigmoid ranked regions
    prior = sigmoid_func(max_rank)

    # Similarities
    # rank factor values for each region
    ranked = rank_regions(data)
    # sigmoid scale ranked values for each factor
    sigmoid = sigmoid_func(ranked)
    # convolution
    conv = np.matmul(sigmoid, sigmoid.T) / sigmoid.shape[1]
    # find top regions and values
    top_regions, top_vals = find_top_regions(conv)
    del conv
    # sqrt of sim values to spread them out followed by sigmoid scaling
    top_vals_sqrt_sig = sigmoid_func(np.sqrt(top_vals), method='full')

    print('Saving...', flush=True)
    # save top 100 of sigmoid euclidean similarities, similar regions, and
    # percentiles for each region
    out_f = h5py.File(out_path, 'w')
    out_f.create_dataset('top_similarity_vals_sqrt_sigmoid',
                         data=top_vals_sqrt_sig)
    out_f.create_dataset('top_similarity_vals', data=top_vals)
    out_f.create_dataset('top_regions', data=top_regions)
    out_f.create_dataset('ranked_data', data=ranked)
    out_f.create_dataset('sigmoid_data', data=sigmoid)
    out_f.create_dataset('priors', data=prior)
    out_f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
