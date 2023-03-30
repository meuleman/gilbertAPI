import numpy as np
import sys
import h5py


def main(ranks_path, thresh, out_path):
    # load all ranks/coordinates
    hdf5_f = h5py.File(ranks_path, 'r')
    percents = hdf5_f['percentiles'][:]
    regions = np.char.decode(hdf5_f['coordinates'][:])
    hdf5_f.close()

    # find top <thresh> regions for each factor and combine
    top_inds = np.unique(np.argsort(-percents, axis=0)[:thresh])

    # find corresponding ranks/coordinates
    top_percents = percents[top_inds]
    top_regions = regions[top_inds]

    # save
    np.savez(out_path, percentiles=top_percents, coordinates=top_regions, original_inds=top_inds)

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
