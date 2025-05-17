import numpy as np
import sys
import h5py


def sigmoid_func_pre(val_arr, scaling_factor=5):
    # max rank
    max_val = np.max(val_arr, axis=0)

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


def sigmoid_func_post(val_arr, scaling_factor=5):

    # max rank
    max_val = np.max(val_arr)

    # scaling between [-scaling_factor, scaling_factor]
    x = (val_arr / (max_val / (scaling_factor * 2))) - scaling_factor
    sig_arr = 1 / (1 + np.exp(-x))
    return sig_arr


def main(base_path, out_path):
    # load the data
    base_hdf5 = h5py.File(base_path, 'r')
    meanmax_arr = base_hdf5['base_arr'][:]
    base_hdf5.close()

    mean_arr = meanmax_arr[:, :18]
    mean_arr = np.where(mean_arr < 0, 0, mean_arr)

    # scale mean scores between 0 and 1
    mean_sig_arr = sigmoid_func_pre(mean_arr)
    # convolution with mean
    mean_conv_mean = np.matmul(mean_sig_arr,
                               mean_sig_arr.T) / len(mean_sig_arr[0])
    # find top 1000 regions
    top_regions, top_vals = find_top_regions(mean_conv_mean)
    del mean_conv_mean
    # sigmoid scale the top 1000 values between 0 and 1
    top_sqrt_vals = np.sqrt(top_vals)
    top_sqrtsig_vals = sigmoid_func_post(top_sqrt_vals)

    print('Saving...', flush=True)
    # save top 100 of sigmoid euclidean similarities, similar regions, and
    # percentiles for each region
    out_f = h5py.File(out_path, 'w')
    out_f.create_dataset('top_similarity_vals_sqrt_sigmoid',
                         data=top_sqrtsig_vals)
    out_f.create_dataset('top_similarity_vals_sqrt', data=top_sqrt_vals)
    out_f.create_dataset('top_similarity_vals', data=top_vals)
    out_f.create_dataset('top_regions', data=top_regions)
    out_f.create_dataset('sigmoid_arr', data=mean_sig_arr)
    out_f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
