import numpy as np
import sys
import h5py


def main(ranks_path, thresh, out_path):
    # load all ranks/coordinates
    hdf5_f = h5py.File(ranks_path, 'r')
    percents = hdf5_f['percentiles'][:]
    regions = np.char.decode(hdf5_f['coordinates'][:])
    hdf5_f.close()

    # find any region index where at least 1 factor rank >= thresh
    inds_set = set()
    for i in range(percents.shape[1]):
        inds = np.where(percents[:, i] >= thresh)[0]
        inds_set.update(inds)

    # find corresponding ranks/coordinates
    top_inds = np.sort(list(inds_set))
    top_percents = percents[top_inds]
    top_regions = regions[top_inds]

    # save
    np.savez(out_path, percentiles=top_percents, coordinates=top_regions, original_inds=top_inds)

if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]), sys.argv[3])
