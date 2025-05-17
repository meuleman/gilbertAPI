import numpy as np
import sys
import tabix
import json
from simsearch_functions import simsearch


def main(DHS_per_f, epi_per_f, gz_file_path,
         regions_path, genomic_location, target_factors):

    # load roi
    gz_tb = tabix.open(gz_file_path)
    roi = np.array(list(gz_tb.querys(genomic_location)))[:, 3].astype(int)[0]

    # find most similar region inds
    (top_regions, roi_sim, roi_conv,
     roi_per, factor_order, metric_order,
     colors) = simsearch(
         DHS_per_f, epi_per_f, roi, target_factors
     )

    # genomic coordinates for each region
    f = open(regions_path, 'r')
    regions = []
    for line in f.readlines():
        line = line.strip().split('\t')
        regions.append(line[0] + ':' + line[1] + '-' + line[2])
    f.close()

    top_region_coords = np.array(regions)[top_regions]

    # for the roi, return results as json array
    # each element is an object of scores for each metric
    return_arr = []
    for rank in range(101):
        region_dict = {
            'rank': rank,
            'coordinates': top_region_coords[rank],
            'region_ind': int(top_regions[rank]),
            'score': float(roi_sim[rank]),
            'factor_conv': list(roi_conv[rank]),
            'percentiles': list(roi_per[rank]),
            'factors': list(factor_order),
            'metrics': list(metric_order),
            'factor_colors': list(colors)
        }
        return_arr.append(region_dict)

    json_return = json.dumps(return_arr)
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
         np.int_(sys.argv[6].strip('[').strip(']').split(',')))
