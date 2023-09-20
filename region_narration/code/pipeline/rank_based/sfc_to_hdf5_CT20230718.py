import numpy as np
import sys
import h5py


def main(mapped_scores_filename, out_path, dataset):
    # load score data
    if dataset == 'cs':
        num_factors = 18
        convert2paper_order = np.arange(18)
    elif dataset == 'dhs':
        num_factors = 16
        convert2paper_order = np.array([8, 3, 5, 1, 4, 15, 9, 6, 2,
                                        12, 13, 11, 0, 10, 7, 14])
    factors = np.loadtxt(mapped_scores_filename, max_rows=1, usecols=np.arange(3, 3 + num_factors), dtype=str, delimiter='\t')[convert2paper_order]
    scores = np.loadtxt(mapped_scores_filename, skiprows=1, usecols=np.arange(3, 3 + num_factors), dtype=float, delimiter='\t')[:, convert2paper_order]

    # set regions with no scores to 0
    scores = np.where(np.isnan(scores), 0, scores)

    # save state means
    ooe_f = h5py.File(out_path, 'w')
    ooe_f.create_dataset('sfc', data=scores)
    ooe_f.create_dataset('factors', data=np.char.encode(factors))
    ooe_f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
