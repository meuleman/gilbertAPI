import numpy as np
import sys
import h5py
import re


def main(mapped_counts_filename, out_path):
    # load count data
    data = []
    f = open(mapped_counts_filename, 'r')
    for line in f.readlines()[1:]:  # not header
        data.append(re.split('[\t|]', line.strip()))
    f.close()
    data = np.array(data)

    # seperating counts, means, and scores
    counts = data[:, 3::3].astype(int)
    means = data[:, 4::3]
    means = np.where(means == 'NAN', np.nan, means).astype('float32')

    # find mean and std for each factor
    factor_means = np.array([])
    factor_std = np.array([])
    for factor in range(len(data[0, 5::3])):
        factor_data = data[:, 5 + (factor * 3)][
                          data[:, 5 + (factor * 3)] != ''
                      ]
        factor_scores = np.array(sum(list(np.char.rsplit(factor_data, ';')),
                                     [])).astype('float32')
        factor_means = np.append(factor_means, np.mean(factor_scores))
        factor_std = np.append(factor_std, np.std(factor_scores))

    # z score
    z = ((means - factor_means) / factor_std * np.sqrt(counts))
    # one sided
    z = np.where(((z < 0) | (np.isnan(z))), 0, z)

    # save count array
    count_f = h5py.File(out_path, 'w')
    count_f.create_dataset('base_arr', data=z)
    count_f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
