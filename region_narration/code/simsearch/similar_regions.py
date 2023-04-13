#!/usr/bin/env python

import numpy as np
import sys
import tabix
import json
import h5py
from simsearch_functions import simsearch, roiURL, load_per


def load_roi(genomic_location, gz_file_path):
    global roi

    # determine roi from hdf5 scale file
    query_start_end = genomic_location.split(':')[1].split('-')
    query_start = int(query_start_end[0])
    query_end = int(query_start_end[1])

    gz_tb = tabix.open(gz_file_path)
    possible_coords = np.array(list(gz_tb.querys(genomic_location)))

    start_dist = np.abs(possible_coords[:, 1].astype(int) - query_start)
    end_dist = np.abs(possible_coords[:, 2].astype(int) - query_end)

    query_region = possible_coords[np.argmin(start_dist + end_dist)]

    roi = int(query_region[3])


def load_ranks(query_path, search_space_path):
    global roi_percents, percents, roi_coords, coords

    # load roi ranks
    hdf5_f = h5py.File(query_path, 'r')
    roi_percents = hdf5_f['percentiles'][roi]
    roi_coords = np.char.decode(hdf5_f['coordinates'][roi])
    hdf5_f.close()

    # load search space ranks
    search_space = np.load(search_space_path)
    percents = search_space['percentiles']
    coords = search_space['coordinates']


def setup_included_factors(included_factors_str, query_factor_thresh):
    global included_factors
    included_factors = np.int_(
        included_factors_str.strip('[').strip(']').split(',')
    )
    query_factors_above_thresh = np.where(roi_percents > query_factor_thresh)[0]
    included_factors = included_factors[
        np.isin(included_factors, query_factors_above_thresh)
    ]


def find_most_similar_region_indices(num_similar_regions):
    global top_regions, top_region_coords, top_regions_sim, top_regions_per, query_factor_order, target_factors, pval, Z, query_sim_prod_ranks
    (top_regions, top_region_coords, top_regions_sim,
     top_regions_per, query_factor_order,
     target_factors, pval,
     Z, query_sim_prod_ranks) = simsearch(roi_percents, percents, roi, roi_coords, coords, included_factors, num_similar_regions)


def setup_roiURLs(get_roiURL, target_factors):
    global roiURLs
    roiURLs = []
    for dl in range(1, len(included_factors) + 1):
        if get_roiURL == 'True':
            # roiURLs.append('roiURL not available')
            roiURLs.append(roiURL(roi, top_region_coords[dl-1],
                                  dl, '_'.join(target_factors.astype(str))))
        else:
            roiURLs.append('roiURL not available')


def reorder_by_query_region():
    global return_top_regions_per, return_pval, return_Z
    return_top_regions_per = top_regions_per[:, :, query_factor_order]
    return_pval = pval[:, query_factor_order]
    return_Z = Z[:, query_factor_order]


def generate_result_array(num_similar_regions):
    global return_arr_full
    return_arr_full = []
    for dl in range(len(included_factors)):
        # for the roi, return results as json array
        # each element is an object of scores for each metric
        return_arr = []
        for rank in range(num_similar_regions + 1):
            coord = top_region_coords[dl, rank]
            region_dict = {
                'rank': rank,
                'coordinates': coord[0] + ':' + coord[1] + '-' + coord[2],
                'region_ind': int(top_regions[dl, rank]),
                'score': float(top_regions_sim[dl, rank]),
                'percentiles': return_top_regions_per[dl, rank].tolist(),
                'prod_rank': query_sim_prod_ranks[dl, rank]
            }
            if rank == 0:
                region_dict['region_ind'] = 'NA'
                region_dict['score'] = 'NA'
                region_dict['query_factor_order'] = query_factor_order.tolist()
                region_dict['target_factors'] = target_factors.tolist()
                region_dict['top_region_factor_Z_score'] = return_Z[dl].tolist()
                region_dict['top_region_factor_pval'] = return_pval[dl].tolist()
                region_dict['roiURL'] = roiURLs[dl]
                region_dict['detail_level'] = dl + 1
                # region_dict['query_prod_rank'] = query_sim_prod_ranks[dl]
            return_arr.append(region_dict)

        return_arr_full.append(return_arr)


def main(query_path, search_space_path, gz_file_path, genomic_location,
         included_factors_str, get_roiURL, num_similar_regions,
         query_factor_thresh):

    # load ROI
    load_roi(genomic_location, gz_file_path)

    # load region percentiles
    load_ranks(query_path, search_space_path)

    # factors the user cares about
    setup_included_factors(included_factors_str, query_factor_thresh)

    # find most similar region inds
    find_most_similar_region_indices(num_similar_regions)

    # dual track roiURL
    setup_roiURLs(get_roiURL, target_factors)

    # reorder based on query region
    reorder_by_query_region()

    # generate result array
    generate_result_array(num_similar_regions)

    json_return = json.dumps(return_arr_full)
    sys.stdout.write(json_return)

    # save output for dev mode
    # f = open('/Users/ctrader/Documents/AG/NarrativeGenomics/' +
    #          'region_narration/code/pipeline/region_similarity_react/' +
    #          'region-narration-app/src/' +
    #          'devmode_top_similar_regions_results_CT20221129.json', 'w')
    # f.write(json_return)
    # f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
         sys.argv[6], int(sys.argv[7]), float(sys.argv[8]))
