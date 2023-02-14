import numpy as np
import sys
import h5py


def load_ranks(ranks_path):
    r = h5py.File(ranks_path, 'r')
    ranks = r['ranked_data'][:]
    max_ranks = r['max_ranks'][:]
    r.close()
    return ranks, max_ranks


def main(DHS_ranks_path, epi_ranks_path, DHS_mapped_path, out_path):
    # load the data
    print('Loading data...', flush=True)

    DHSalpha2paper_order = np.array([8, 3, 5, 1, 4, 15, 9, 6, 2, 12, 13,
                                     11, 0, 10, 7, 14])

    DHS_ranks, DHS_max_ranks = load_ranks(DHS_ranks_path)
    # reorder DHS factors to paper order
    DHS_ranks = DHS_ranks[:, DHSalpha2paper_order]
    DHS_max_ranks = DHS_max_ranks[DHSalpha2paper_order]

    epi_ranks, epi_max_ranks = load_ranks(epi_ranks_path)

    epi_factor_order = np.array(['TssA', 'TssFlnk', 'TssFlnkU',
                                 'TssFlnkD', 'Tx', 'TxWk', 'EnhG1',
                                 'EnhG2', 'EnhA1', 'EnhA2', 'EnhWk',
                                 'ZNF_Rpts', 'Het', 'TssBiv', 'EnhBiv',
                                 'ReprPC', 'ReprPCWk', 'Quies'])

    f = open(DHS_mapped_path, 'r')
    DHS_factor_order = np.array(
        f.readline().strip().split('\t')
    )[3:][DHSalpha2paper_order]
    f.close()

    # join metrics
    joined_ranks = np.hstack((DHS_ranks, epi_ranks))
    joined_max_ranks = np.append(DHS_max_ranks, epi_max_ranks)
    joined_factor_order = np.append(DHS_factor_order, epi_factor_order)
    joined_metric_order = np.append(
        np.array(['DHS signal'] * len(DHS_factor_order)),
        np.array(['chromatin states'] * len(epi_factor_order))
    )

    percents = joined_ranks / joined_max_ranks

    # for each region, convolve against all other regions
    print('Convolution...', flush=True)
    conv_top = []
    sim_top = []
    regions_top = []
    for i in range(len(percents)):
        # convolve
        i_percents = percents[i]
        i_conv_all = percents * i_percents
        # fusion similarity function
        i_sim_all = 1 - np.prod(1 - i_conv_all, axis=1)
        # find top scoring regions
        i_top_regions = np.argsort(-i_sim_all)[:101]
        i_top_regions = i_top_regions[i_top_regions != i][:100]
        i_conv_top = i_conv_all[i_top_regions]
        i_sim_top = i_sim_all[i_top_regions]
        conv_top.append(i_conv_top)
        sim_top.append(i_sim_top)
        regions_top.append(i_top_regions)

    conv_top = np.array(conv_top)
    sim_top = np.array(sim_top)
    regions_top = np.array(regions_top)

    print('Saving...', flush=True)
    # save top 100 similarity scores, similar regions, and
    # convolution scores for each region
    out_f = h5py.File(out_path, 'w')
    out_f.create_dataset('convolution', data=conv_top)
    out_f.create_dataset('similarity', data=sim_top)
    out_f.create_dataset('regions', data=regions_top)
    out_f.create_dataset('factor_order', data=joined_factor_order)
    out_f.create_dataset('metric_order', data=joined_metric_order)
    out_f.create_dataset('ranks', data=joined_ranks)
    out_f.create_dataset('max_ranks', data=joined_max_ranks)
    out_f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
