import numpy as np
from scipy.stats import rankdata


"""
Takes in a region x factor array of scores. Regions are ranked
relative to each other and converted to percentile values
between 0 and 1.
"""


def region_percentile(arr):
    # rank each region
    ranked_regions = rankdata(arr, axis=0, method='min') - 1

    # # find the number of regions with 0 values for each component
    # comp_num_zero = len(arr) - np.count_nonzero(arr, axis=0)

    # # adjusting/resetting ranks
    # ranked_regions = ranked_regions - comp_num_zero
    # ranked_regions = np.where(ranked_regions < 0, 0, ranked_regions)

    # percentile = ranked_regions / np.max(ranked_regions, axis=0)
    percentile = ranked_regions / (ranked_regions.shape[0] - 1)

    return percentile
