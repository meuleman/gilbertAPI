import numpy as np
import sys
import h5py
from rankbased_preprocessing_functions import region_percentile


def main(sum_path, out_path):
    # convert alphabetical ordering to official ordering
    epialpha2paper_order = np.array([10, 12, 14, 13, 15, 16, 3, 4, 0, 1,
                                     5, 17, 6, 11, 2, 8, 9, 7])

    # load factor sums (in official order)
    hdf5_f = h5py.File(sum_path, 'r')
    sums = hdf5_f['state_sums'][:][:, epialpha2paper_order]
    factors = np.char.decode(hdf5_f['factors'][:])[epialpha2paper_order]
    hdf5_f.close()
    sums = np.where(sums < 0, 0, sums)

    epi_color = np.array(['#ff0000', '#ff4500', '#ff4500', '#ff4500',
                          '#008000', '#006400', '#c2e105', '#c2e105',
                          '#ffc34d', '#ffc34d', '#ffff00', '#66cdaa',
                          '#8a91d0', '#cd5c5c', '#bdb76b', '#808080',
                          '#c0c0c0', '#ffffff'])

    percentile = region_percentile(sums)

    # save
    ooe_f = h5py.File(out_path, 'w')
    ooe_f.create_dataset('percentile', data=percentile)
    ooe_f.create_dataset('factors', data=np.char.encode(factors))
    ooe_f.create_dataset('colors', data=np.char.encode(epi_color))
    ooe_f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
