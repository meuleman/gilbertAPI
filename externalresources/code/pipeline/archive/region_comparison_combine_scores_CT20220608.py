import numpy as np
import sys
import os
import time


def combine_scores(DHS_count_rank_vals, DHS_count_similar_regions, TF_count_rank_vals, TF_count_similar_regions, epilogos_score_rank_vals, epilogos_score_similar_regions):
    #max ind for combined score array
    max_ind = np.max([DHS_count_similar_regions, TF_count_similar_regions, epilogos_score_similar_regions])
    #DHS array
    DHS_arr = np.zeros(max_ind + 1)
    DHS_arr[DHS_count_similar_regions] = top_100_rank_vals
    return

def main(region_of_interest_path, DHS_npz_path, TF_npz_path, epilogos_npz_path):
    t = time.time()
    #find the region of interest
    f = open(region_of_interest_path, 'r')
    roi = int(f.readline())

    #####Load rank and ind arrays (np.load('')['this array'][roi]
    #DHS counts
    DHS_count_rank_vals = np.load(DHS_npz_path)['top_1000_rank_vals'][roi]
    DHS_count_similar_regions = np.load(DHS_npz_path)['top_1000_similar_regions'][roi]

    #TF motif counts
    TF_count_rank_vals = np.load(TF_npz_path)['top_1000_rank_vals'][roi]
    TF_count_similar_regions = np.load(TF_npz_path)['top_1000_similar_regions'][roi]

    #epilogos score sum and max
    epilogos_score_rank_vals = np.load(epilogos_npz_path)['top_1000_rank_vals'][roi]
    epilogos_score_similar_regions = np.load(epilogos_npz_path)['top_1000_similar_regions'][roi]

    combine_scores(DHS_count_rank_vals, DHS_count_similar_regions, TF_count_rank_vals, TF_count_similar_regions, epilogos_score_rank_vals, epilogos_score_similar_regions)

    print(time.time()-t)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
