import numpy as np
from functools import partial
from multiprocessing import Pool
from contextlib import closing
import sys
import h5py
from collections import Counter
import re
from scipy.stats import fisher_exact

def load_hdf5_files_full(hdf5_path, max_ind, roi):
    #take in a file path and load scaled scores and similar regions
    hdf5_file = h5py.File(hdf5_path, 'r')
    scaled_vals = hdf5_file['top_1000_scaled_vals'][roi]
    similar_regions = hdf5_file['top_1000_similar_regions'][roi]
    hdf5_file.close()
    arr = np.zeros(max_ind + 1)
    arr[similar_regions] = scaled_vals
    return arr

def find_scores(i):
    DHS_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220811.JQ.DHScosine.hdf5'
    TF_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220811.JQ.TFcosine.hdf5'
    epilogos_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220811.JQ.epilogoscosine.hdf5'
    coex_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220811.JQ.coexpression_human_transferred_maxSTRINGsimilarity.hdf5'
    cooccur_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220811.JQ.cooccurenceSTRINGsimilarity.hdf5'
    database_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220811.JQ.database_human_transferred_maxSTRINGsimilarity.hdf5'
    exper_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220811.JQ.experiments_human_transferred_maxSTRINGsimilarity.hdf5'
    fusion_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220811.JQ.fusionSTRINGsimilarity.hdf5'
    homology_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220811.JQ.homologySTRINGsimilarity.hdf5'
    neigh_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220811.JQ.neighborhood_human_transferred_maxSTRINGsimilarity.hdf5'
    textmin_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220811.JQ.textmining_human_transferred_maxSTRINGsimilarity.hdf5'

    max_ind=100001
    DHS_arr = load_hdf5_files_full(DHS_hdf5_path, max_ind, i)
    TF_arr = load_hdf5_files_full(TF_hdf5_path, max_ind, i)
    epilogos_arr = load_hdf5_files_full(epilogos_hdf5_path, max_ind, i)
    coex_arr = load_hdf5_files_full(coex_hdf5_path, max_ind, i)
    cooccur_arr = load_hdf5_files_full(cooccur_hdf5_path, max_ind, i)
    database_arr = load_hdf5_files_full(database_hdf5_path, max_ind, i)
    exper_arr = load_hdf5_files_full(exper_hdf5_path, max_ind, i)
    fusion_arr = load_hdf5_files_full(fusion_hdf5_path, max_ind, i)
    homology_arr = load_hdf5_files_full(homology_hdf5_path, max_ind, i)
    neigh_arr = load_hdf5_files_full(neigh_hdf5_path, max_ind, i)
    textmin_arr = load_hdf5_files_full(textmin_hdf5_path, max_ind, i)

    #scoring all metrics
    score_arr_i = 1 - ((1 - epilogos_arr) * (1 - DHS_arr) * (1 - TF_arr) * (1 - coex_arr) * (1 - cooccur_arr) * (1 - database_arr) * (1 - exper_arr) *
                     (1 - fusion_arr) * (1 - homology_arr) * (1 - neigh_arr) * (1 - textmin_arr))
    top_100_inds = np.argsort(-score_arr_i)[:100]
    scores = np.vstack((score_arr_i[top_100_inds], DHS_arr[top_100_inds], TF_arr[top_100_inds], 
                        epilogos_arr[top_100_inds], coex_arr[top_100_inds], cooccur_arr[top_100_inds], 
                        database_arr[top_100_inds], exper_arr[top_100_inds], fusion_arr[top_100_inds], 
                        homology_arr[top_100_inds], neigh_arr[top_100_inds], textmin_arr[top_100_inds]))
    return scores

def main(num_cpus):
    #split on cluster and find significant terms
    print('Scores for all regions...')
    func = partial(find_scores)
    with closing(Pool(num_cpus)) as pool:
        scores = pool.map(func, np.arange(98629))
    pool.join()
    scores = np.array(scores)

    np.save('/home/ctrader/region_similarity/data/regions_25kb_CT20220812.JQ.scores_for_all_regions.fusion_and_11_metrics.npy', scores)

if __name__ == "__main__":
    main(int(sys.argv[1]))
