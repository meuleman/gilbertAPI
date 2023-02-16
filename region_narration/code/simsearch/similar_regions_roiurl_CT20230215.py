import numpy as np
import sys
import tabix
from simsearch_functions import simsearch, roiURL


def main(DHS_per_f, epi_per_f, gz_file_path,
         regions_path, genomic_location, target_factors):

    # load roi
    gz_tb = tabix.open(gz_file_path)
    roi = np.array(list(gz_tb.querys(genomic_location)))[:, 3].astype(int)[0]

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
        regions.append(line)
    f.close()

    top_region_coords = np.array(regions)[top_regions]

    roiurl_in = np.hstack((top_region_coords, roi_sim.reshape(-1, 1)))
    print(roiURL(roi, roiurl_in))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
         np.int_(sys.argv[6].strip('[').strip(']').split(',')))
