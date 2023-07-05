import numpy as np
import sys
import h5py


def main(ranks_path, thresh, out_path):
    # load all ranks/coordinates
    hdf5_f = h5py.File(ranks_path, 'r')
    ranks = hdf5_f['ranks'][:]
    regions = np.char.decode(hdf5_f['coordinates'][:])
    hdf5_f.close()

    # find top <thresh> regions for each factor (with no overlap) and combine
    overall_top_inds = []
    for factor in range(ranks.shape[1]):
        factor_top_inds = []
        factor_top_regions = (-1 * np.ones((thresh, 3), dtype=int)).astype(str)

        factor_sorted_inds = np.argsort(-ranks[:, factor])
        factor_sorted_regions = regions[factor_sorted_inds]

        regions_added = 0
        for ind, coords in zip(factor_sorted_inds, factor_sorted_regions):
            num_overlaps = np.count_nonzero(
                (factor_top_regions[:, 0].astype(str) == coords[0]) &
                (
                    (
                        (factor_top_regions[:, 1].astype(int) > int(coords[1])) &
                        (factor_top_regions[:, 1].astype(int) < int(coords[2]))
                    ) |
                    (
                        (factor_top_regions[:, 2].astype(int) > int(coords[1])) &
                        (factor_top_regions[:, 2].astype(int) < int(coords[2]))
                    ) |
                    (
                        (factor_top_regions[:, 1].astype(int) == int(coords[1])) &
                        (factor_top_regions[:, 2].astype(int) == int(coords[2]))
                    )
                )
            )
            if num_overlaps == 0:
                factor_top_inds.append(ind)
                factor_top_regions[regions_added] = coords
                regions_added += 1
            if len(factor_top_inds) == thresh:
                break
        overall_top_inds.append(factor_top_inds)
    overall_top_inds = np.array(overall_top_inds) 

    top_inds = np.unique(overall_top_inds)

    # find corresponding ranks/coordinates
    top_ranks = ranks[top_inds]
    top_regions = regions[top_inds]

    # find the factors each ind was selected for
    factors_for_top_inds = []
    for ind in top_inds:
        factors_for_ind = np.where(overall_top_inds == ind)[0]
        factors_for_top_inds.append(factors_for_ind)
    factors_for_top_inds = np.array(factors_for_top_inds, dtype=object)

    # save
    np.savez(out_path, ranks=top_ranks, coordinates=top_regions, original_inds=top_inds, factors=factors_for_top_inds)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
