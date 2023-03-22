#!/usr/bin/env python

import numpy as np
import sys
import tabix
import json
from simsearch_functions import simsearch, roiURL, load_per

def load_roi(genomic_location, gz_file_path):
    global roi
    if str.isalpha(str(genomic_location)[0]):
        # load roi
        gz_tb = tabix.open(gz_file_path)
        roi = np.array(
            list(gz_tb.querys(genomic_location))
        )[:, 3].astype(int)[0]
    else:
        roi = int(genomic_location)


def load_ranks(per_f):
    global percents
    percents = np.load(per_f)


def setup_included_factors(included_factors_str):
    global included_factors
    included_factors = np.int_(
        included_factors_str.strip('[').strip(']').split(',')
    )


def read_in_genomic_coordinates(regions_path):
    global regions_sep
    regions_sep = np.load(regions_path)


def find_most_similar_region_indices(num_similar_regions):
    global top_regions, roi_sim, roi_per, query_factor_order, target_factors, pval, Z
    (top_regions, roi_sim,
     roi_per, query_factor_order,
     target_factors, pval,
     Z) = simsearch(percents, roi, included_factors, num_similar_regions)


def setup_coordinates(top_regions):
    global top_region_coords
    top_region_coords = regions_sep[top_regions]


def setup_roiURLs(get_roiURL, top_regions, target_factors):
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
    global return_roi_per, return_pval, return_Z
    return_roi_per = roi_per[:, :, query_factor_order]
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
                'score': float(roi_sim[dl, rank]),
                'percentiles': return_roi_per[dl, rank].tolist(),
            }
            if rank == 0:
                region_dict['query_factor_order'] = query_factor_order.tolist()
                region_dict['target_factors'] = target_factors.tolist()
                region_dict['top_region_factor_Z_score'] = return_Z[dl].tolist()
                region_dict['top_region_factor_pval'] = return_pval[dl].tolist()
                region_dict['roiURL'] = roiURLs[dl]
                region_dict['detail_level'] = dl + 1
            return_arr.append(region_dict)

        return_arr_full.append(return_arr)


def main(per_f, gz_file_path, regions_path, genomic_location,
         included_factors_str, get_roiURL, num_similar_regions):

    # load ROI
    load_roi(genomic_location, gz_file_path)

    # load region percentiles
    load_ranks(per_f)

    # factors the user cares about
    setup_included_factors(included_factors_str)

    # read in genomic coordinates for each region
    read_in_genomic_coordinates(regions_path)

    # find most similar region inds
    find_most_similar_region_indices(num_similar_regions)

    # get coordinates
    setup_coordinates(top_regions)

    # dual track roiURL
    setup_roiURLs(get_roiURL, top_regions, target_factors)

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
         sys.argv[6], int(sys.argv[7]))
