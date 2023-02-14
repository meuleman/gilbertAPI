import numpy as np
import sys
import h5py
from rankbased_preprocessing_functions import region_percentile


def main(NMF_path, out_path):
    # convert alphabetical ordering to official ordering
    DHSalpha2paper_order = np.array([8, 3, 5, 1, 4, 15, 9, 6, 2,
                                     12, 13, 11, 0, 10, 7, 14])

    # load NMF sums (in official order)
    hdf5_f = h5py.File(NMF_path, 'r')
    factors = np.char.decode(hdf5_f['factors'][:])[DHSalpha2paper_order]
    NMFsum = hdf5_f['factor_sums'][:][:, DHSalpha2paper_order]
    hdf5_f.close()

    NMFsum_total = np.sum(NMFsum, axis=1)
    NMFenrich = NMFsum / NMFsum_total.reshape(-1, 1)
    NMFenrich = np.where(np.isnan(NMFenrich), 0, NMFenrich)

    DHS_color = np.array(['#FFE500', '#FE8102', '#FF0000', '#07AF00',
                          '#4C7D14', '#414613', '#05C1D9', '#0467FD',
                          '#009588', '#BB2DD4', '#7A00FF', '#4A6876',
                          '#08245B', '#B9461D', '#692108', '#C3C3C3'])

    # enrichment and sum percentiles
    NMFenrich_per = region_percentile(NMFenrich)
    NMFsum_per = region_percentile(NMFsum)
    # co-rank enrichment and sum percentiles
    percentile = (NMFenrich_per + NMFsum_per) / 2

    # save
    ooe_f = h5py.File(out_path, 'w')
    ooe_f.create_dataset('percentile', data=percentile)
    ooe_f.create_dataset('factors', data=np.char.encode(factors))
    ooe_f.create_dataset('colors', data=np.char.encode(DHS_color))
    ooe_f.create_dataset('NMFenrichment_percentile', data=NMFenrich_per)
    ooe_f.create_dataset('NMFsum_percentile', data=NMFsum_per)
    ooe_f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
