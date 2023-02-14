import numpy as np
import sys
import h5py
import re


def main(mapped_counts_filename, out_path):
    # scores from all regions
    counts = []
    means = []
    scores = np.array([])

    f = open(mapped_counts_filename, 'r')

    # array to track scores for each factor
    factors = f.readline().strip().split('\t')[3:]

    # for each region
    for line in f.readlines():  # not header
        line = re.split('[\t|]', line.strip())

        # inidvidual scores for each region
        line_scores_str = line[5::3]
        line_scores_str = np.char.rsplit(line_scores_str, ';')
        line_scores = np.array([])
        for factor in range(len(factors)):
            # test to see if first element is not ''
            if line_scores_str[factor][0]:
                line_scores = np.append(line_scores, np.float_(line_scores_str[factor]))
        # number of factors for region
        counts.append(len(line_scores))
        # average score for region
        if len(line_scores) > 0:
            means.append(np.mean(line_scores))
        else:
            means.append(0)
        # scores from all regions
        scores = np.append(scores, line_scores)
    f.close()

    counts = np.array(counts)
    means = np.array(means)

    # find mean and std for each factor
    gw_mean = np.mean(scores)
    gw_std = np.std(scores)

    # z score
    z = ((means - gw_mean) / gw_std * np.sqrt(counts))
    # one sided
    z = np.where(z < 0, 0, z)

    # save count array
    count_f = h5py.File(out_path, 'w')
    count_f.create_dataset('factor_blind_z', data=z)
    count_f.create_dataset('counts', data=counts)
    count_f.create_dataset('means', data=means)
    count_f.create_dataset('gw_mean', data=gw_mean)
    count_f.create_dataset('gw_std', data=gw_std)
    count_f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
