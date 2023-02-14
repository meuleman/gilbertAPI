import numpy as np
import sys
import h5py
import re


def main(mapped_counts_filename, out_path):
    # load count data
    counts = []
    means = []

    f = open(mapped_counts_filename, 'r')

    # array to track scores for each factor
    factors = f.readline().strip().split('\t')[3:]
    factor_scores = []
    for i in range(len(factors)):
        factor_scores.append([])

    # for each region
    for line in f.readlines():  # not header
        line = re.split('[\t|]', line.strip())
        # number of each motif for region
        counts.append(line[3::3])
        # average score for each motif for region
        means.append(line[4::3])

        # inidvidual scores for each region
        line_scores = line[5::3]
        line_scores = np.char.rsplit(line_scores, ';')
        for factor in range(len(factors)):
            # test to see if first element is ''
            if line_scores[factor][0]:
                factor_scores[
                    factor
                ].extend(list(np.float_(line_scores[factor])))
    f.close()

    counts = np.int_(counts)
    means = np.array(means)
    means = np.where(means == 'NAN', np.nan, means).astype('float32')

    # find mean and std for each factor
    factor_means = np.array([])
    factor_std = np.array([])
    for factor in range(len(factors)):
        factor_means = np.append(factor_means, np.mean(factor_scores[factor]))
        factor_std = np.append(factor_std, np.std(factor_scores[factor]))

    # z score
    z = ((means - factor_means) / factor_std * np.sqrt(counts))
    # one sided
    z = np.where(((z < 0) | (np.isnan(z))), 0, z)

    # save count array
    count_f = h5py.File(out_path, 'w')
    count_f.create_dataset('z', data=z)
    count_f.create_dataset('means', data=means)
    count_f.create_dataset('counts', data=counts)
    count_f.create_dataset('factor_means', data=factor_means)
    count_f.create_dataset('factor_std', data=factor_std)
    count_f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
