import numpy as np
import sys
import h5py


def main(mapped_scores_filename, out_path):
    # load score data
    scores = []
    f = open(mapped_scores_filename, 'r')
    factors = f.readline().strip().split('\t')[3:]
    for line in f.readlines():
        scores.append(line.strip().split('\t')[3:])
    f.close()
    scores = np.float_(scores)

    # save state means
    ooe_f = h5py.File(out_path, 'w')
    ooe_f.create_dataset('state_means', data=scores)
    ooe_f.create_dataset('factors', data=np.char.encode(factors))
    ooe_f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
