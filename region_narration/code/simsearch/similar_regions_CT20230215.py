import numpy as np
import os
# from datetime import date
import sys
import tabix
import json
from simsearch_functions import simsearch, roiURL

def main(DHS_per_f, epi_per_f, gz_file_path,
         regions_path, genomic_location, target_factors):

    # load roi
    gz_tb = tabix.open(gz_file_path)
    roi = np.array(list(gz_tb.querys(genomic_location)))[:, 3].astype(int)

    # find most similar region inds
    (top_regions, roi_sim, roi_conv,
     factor_order, metric_order, colors) = simsearch(
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

    roiurl_in = np.hstack((top_region_coords, roi_sim.reshape(-1, 1)))
    print(roiURL(roi, roiurl_in))


    #sort_order = np.argsort(top_region_inds[0])
    ## sentences for all regions
    #sentences = make_sentences(DHS_hdf5_path, DHS_mapped_path, TF_hdf5_path,
    #                           TF_mapped_path, epilogos_hdf5_path,
    #                           epilogos_mapped_path, STRING_ranks_hdf5_path,
    #                           roi=top_region_inds[0][sort_order])
    #sentences = sentences[np.argsort(sort_order)]

    ## for the first roi, return results as json array
    ## each element is an object of scores for each metric
    #chromatin_return = []
    #for rank in range(101):
    #    region_dict = {}
    #    for metric, name in zip(range(len(chromatin_metric_order)),
    #                            chromatin_metric_order):
    #        region_dict[name] = chromatin_arr_top_100[metric, 0, rank]
    #    region_dict['rank'] = rank
    #    region_dict['coordinates'] = top_region_coords[0, rank]
    #    region_dict['region_ind'] = int(top_region_inds[0, rank])
    #    region_dict['sentence'] = sentences[rank]
    #    chromatin_return.append(region_dict)

    #json_return = json.dumps(chromatin_return)
    #sys.stdout.write(json_return)

    ## save output for dev mode
    ## f = open('/Users/ctrader/Documents/AG/NarrativeGenomics/' +
    ##          'region_narration/code/pipeline/region_similarity_react/' +
    ##          'region-narration-app/src/' +
    ##          'devmode_top_similar_regions_results_CT20221129.json', 'w')
    ## f.write(json_return)
    ## f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
         sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10],
         sys.argv[11])
