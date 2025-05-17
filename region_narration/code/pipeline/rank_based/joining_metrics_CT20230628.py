import numpy as np
import sys
import h5py

def load_ranks(ranks_f):
    hdf5_f = h5py.File(ranks_f, 'r')
    ranks = hdf5_f['ranks'][:]
    hdf5_f.close()
    return ranks


def main(dhs_path, epi_path, regions_path, out_path, out_text):
    # load region ranks
    DHS_ranks = load_ranks(dhs_path)
    epi_ranks = load_ranks(epi_path)

    # join percents
    ranks = np.hstack((DHS_ranks, epi_ranks))

    # load region coordinates
    f = open(regions_path, 'r')
    regions = []
    for line in f.readlines():
        line = line.strip().split('\t')
        regions.append(line)
    f.close()
    regions = np.array(regions)

    # save
    ranks_f = h5py.File(out_path, 'w')
    ranks_f.create_dataset('ranks', data=ranks.astype('float32'))
    ranks_f.create_dataset('coordinates', data=np.char.encode(regions))
    ranks_f.close()

    # save txt file for indexing
    np.savetxt(out_text, ranks.astype('float32'), fmt='%s', delimiter='\t')

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
