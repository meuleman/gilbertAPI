import numpy as np
import sys
import h5py


def main(mapped_counts_filename, out_path):
    # load count data
    count_data = []
    f = open(mapped_counts_filename, 'r')
    for line in f:
        count_data.append(line.strip().split('\t'))
    f.close()
    count_data = np.array(count_data)

    # seperating count labels and count values
    factors = count_data[0][3:]
    count_data_vals = count_data[1:, 3:].astype(int)

    # sum counts in region
    region_sum = np.sum(count_data_vals, axis=1)
    # pseudocount for regions with no factors
    region_sum = np.where(region_sum == 0, 1, region_sum)

    # ratio of factor counts / total counts in region
    obs_ratio = (count_data_vals.T / region_sum).T

    # save count array
    count_f = h5py.File(out_path, 'w')
    count_f.create_dataset('counts', data=count_data_vals)
    count_f.create_dataset('factors', data=np.char.encode(factors))
    count_f.create_dataset('enrichment', data=obs_ratio)
    count_f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
