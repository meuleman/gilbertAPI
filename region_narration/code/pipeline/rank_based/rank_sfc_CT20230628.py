import numpy as np
import sys
import h5py
from rankbased_preprocessing_functions import region_percentile


def main(sfc_path, out_path, dataset):
    # convert alphabetical ordering to official ordering
    if dataset == 'cs':
        # convert2paper_order = np.array([10, 12, 14, 13, 15, 16, 3, 4, 0, 1,
        #                                 5, 17, 6, 11, 2, 8, 9, 7])
        convert2paper_order = np.arange(18)

        colors = np.array(['#ff0000', '#ff4500', '#ff4500', '#ff4500',
                           '#008000', '#006400', '#c2e105', '#c2e105',
                           '#ffc34d', '#ffc34d', '#ffff00', '#66cdaa',
                           '#8a91d0', '#cd5c5c', '#bdb76b', '#808080',
                           '#c0c0c0', '#ffffff'])

    elif dataset == 'dhs':
        convert2paper_order = np.array([8, 3, 5, 1, 4, 15, 9, 6, 2,
                                        12, 13, 11, 0, 10, 7, 14])

        colors = np.array(['#FFE500', '#FE8102', '#FF0000', '#07AF00',
                           '#4C7D14', '#414613', '#05C1D9', '#0467FD',
                           '#009588', '#BB2DD4', '#7A00FF', '#4A6876',
                           '#08245B', '#B9461D', '#692108', '#C3C3C3'])

    # load factor sfc scores (in official order)
    hdf5_f = h5py.File(sfc_path, 'r')
    sfc = hdf5_f['sfc'][:][:, convert2paper_order]
    factors = np.char.decode(hdf5_f['factors'][:])[convert2paper_order]
    hdf5_f.close()
    sfc = np.where(sfc < 0, 0, sfc)

    ranks = region_percentile(sfc)

    # save
    ooe_f = h5py.File(out_path, 'w')
    ooe_f.create_dataset('ranks', data=ranks)
    ooe_f.create_dataset('factors', data=np.char.encode(factors))
    ooe_f.create_dataset('colors', data=np.char.encode(colors))
    ooe_f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
