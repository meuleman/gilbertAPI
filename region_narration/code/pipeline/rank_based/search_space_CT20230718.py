import numpy as np
import sys
import h5py


"""
Function to load SFC values from each dataset
"""
def load_sfc(sfc_path):
    hdf5_f = h5py.File(sfc_path, 'r')
    sfc = hdf5_f['sfc'][:]
    hdf5_f.close()
    sfc = np.where(sfc < 0, 0, sfc)
    return sfc


def main(dhs_sfc_path, cs_sfc_path, ranks_path, top_region_thresh, multi_factor_thresh, out_path):
    # load sfc data
    dhs_sfc = load_sfc(dhs_sfc_path)
    cs_sfc = load_sfc(cs_sfc_path)
    sfc = np.hstack((dhs_sfc, cs_sfc))
    del dhs_sfc, cs_sfc

    # calculate the number of factors significant in each segment
    segment_counts = np.count_nonzero(sfc, axis=1)

    # include any segment with at least <multi_factor_thresh> counts
    segments_to_include = np.where(segment_counts >= multi_factor_thresh)[0]

    # segments we may include
    potential_segments_to_include = np.where(
        (segment_counts < multi_factor_thresh) & 
        (segment_counts > 0)
    )[0]

    # sfc cutoff
    sfc_thresh = np.sort(sfc, axis=0)[-top_region_thresh]

    # potential segments greater than sfc threshold
    mask = (
        (sfc[potential_segments_to_include] >= sfc_thresh) &
        (sfc[potential_segments_to_include] > 0)
    )

    # limit potential segments
    potential_segments_to_include = potential_segments_to_include[
        np.where(
            np.count_nonzero(mask, axis=1) > 0
        )[0]
    ]

    del sfc

    # full list of segments to include
    top_inds = np.unique(
        np.append(segments_to_include, potential_segments_to_include)
    )

    # load all ranks/coordinates
    hdf5_f = h5py.File(ranks_path, 'r')
    ranks = hdf5_f['ranks'][:]
    regions = np.char.decode(hdf5_f['coordinates'][:])
    hdf5_f.close()

    # find corresponding ranks/coordinates
    top_ranks = ranks[top_inds]
    top_regions = regions[top_inds]

    # save
    np.savez(out_path, ranks=top_ranks, coordinates=top_regions, original_inds=top_inds)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), sys.argv[6])
