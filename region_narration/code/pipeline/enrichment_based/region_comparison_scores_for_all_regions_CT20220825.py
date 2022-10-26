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
    DHS_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220815.DHScosine.hdf5'
    TF_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220815.TFcosine.hdf5'
    epilogos_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220815.epilogoscosine.hdf5'
    coex_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220815.coexpression_human_transferred_maxSTRINGsimilarity.hdf5'
    cooccur_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220815.cooccurenceSTRINGsimilarity.hdf5'
    database_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220815.database_human_transferred_maxSTRINGsimilarity.hdf5'
    exper_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220815.experiments_human_transferred_maxSTRINGsimilarity.hdf5'
    fusion_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220815.fusionSTRINGsimilarity.hdf5'
    homology_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220815.homologySTRINGsimilarity.hdf5'
    neigh_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220815.neighborhood_human_transferred_maxSTRINGsimilarity.hdf5'
    textmin_hdf5_path = '/home/ctrader/region_similarity/data/regions_25kb_CT20220815.textmining_human_transferred_maxSTRINGsimilarity.hdf5'

    max_ind=num_regions
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
                        homology_arr[top_100_inds], neigh_arr[top_100_inds], textmin_arr[top_100_inds], 
                        top_100_inds))
    return scores

def main(num_cpus):
    #get number of regions
    global num_regions
    filename = '/home/ctrader/region_similarity/data/Altius_25kb_hg38_maxMean_CT20220815.c1t3.bed'
    f = open(filename, 'r')
    num_regions = len(f.readlines())

    #split on cluster and find significant terms
    print('Scores for all regions...', flush=True)
    func = partial(find_scores)
    with closing(Pool(num_cpus)) as pool:
        scores = pool.map(func, np.arange(num_regions))
    pool.join()
    scores = np.array(scores)
    top_100_regions = scores[:, -1, :].astype(int)
    scores = scores[:, :-1, :]
    top_100_regions = np.hstack((np.arange(len(top_100_regions)).reshape(-1, 1), top_100_regions))

    np.savez('/home/ctrader/region_similarity/data/regions_25kb_CT20220825.scores_and_top_recs_for_all_regions.fusion_and_11_metrics.top_recs.npz', 
             scores=scores, roi_and_top_100_regions=top_100_regions)
if __name__ == "__main__":
    main(int(sys.argv[1]))
