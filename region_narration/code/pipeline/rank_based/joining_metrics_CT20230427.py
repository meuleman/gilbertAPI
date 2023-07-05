import numpy as np
import sys
import h5py

def load_per(per_f):
    hdf5_f = h5py.File(per_f, 'r')
    percentile = hdf5_f['percentile'][:]
    hdf5_f.close()
    return percentile


def main(dhs_path, epi_path, regions_path, out_path, out_text):
    # load region percentiles
    DHS_per = load_per(dhs_path)
    epi_per = load_per(epi_path)

    # join percents
    percents = np.hstack((DHS_per, epi_per))

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
    ranks_f.create_dataset('percentiles', data=percents.astype('float32'))
    ranks_f.create_dataset('coordinates', data=np.char.encode(regions))
    ranks_f.close()

    # save txt file for indexing
    np.savetxt(out_text, percents.astype('float32'), fmt='%s', delimiter='\t')

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
