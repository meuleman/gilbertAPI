from flask import Blueprint, request
from common import hilbert_pos_to_order
import numpy as np
import h5py
from scipy.stats import rankdata
import pickle
# import os
import time
import json
from constants import order_independent_filtering_path, CHR_SIZES_UCSC, precalculated_region_set_enrichment_orders as precalculated_orders

# Create a Blueprint
region_set_enrichment_bp = Blueprint('region_set_enrichment_bp', __name__)


"""
Loads slice of coverage arrays for given region
"""
def load_slice_coverage(chrm, start, end):
    with h5py.File(order_independent_filtering_path, 'r') as f:
        coverage = f[chrm][start:end]

    return np.sum(coverage, axis=0)


"""
Loads array of factors (in order of the coverage columns)
"""
def load_all_factors():
    with h5py.File(order_independent_filtering_path, 'r') as f:
        factors = np.char.decode(f['factors'][:])
    return factors


"""
Loads all coverage arrays for a given region
"""
def load_gw_coverage():
    with h5py.File(order_independent_filtering_path, 'r') as f:
        coverage = f['coverage'][:]
    return coverage


"""
Load coverage counts for each factor for segment at order. Each row contains
the counts for a full segment at the given order.
"""
def load_precalculated_factor_coverage(chrm, order, i):
    with h5py.File(order_independent_filtering_path, 'r') as f:
        coverage = f[f"{chrm}.{order}.counts"][i]
    return coverage


"""
Calculates the enrichment of bp coverage for each factor within a region set. Returns 
the log2 of the observed over expected coverage for each factor.
"""
def get_region_set_factor_enrichments(regions, allow_precalculated=True, membership=False):
    start = time.time()
    # load the genome-wide coverage of each factor
    gw_factor_coverage = load_gw_coverage()

    # for each region in region set, determine coverage of each factor
    if membership:
        region_membership = np.zeros((len(gw_factor_coverage), len(regions)), dtype=bool)
    else:
        region_coverage = 0
        factor_coverage = np.zeros(len(gw_factor_coverage))
        segment_bin_count = np.zeros(len(gw_factor_coverage), dtype=int)

    for i, r in enumerate(regions):
        r_chrm = r['chromosome']
        r_order = r["order"]
        r_i = r["i"]
        r_start = hilbert_pos_to_order(r_i, r_order, 14)
        # cap to end of chromosome for proper normalization
        r_end = min(hilbert_pos_to_order(r_i + 1, r_order, 14), CHR_SIZES_UCSC[r_chrm])

        # collect coverage for each factor within region
        if r_order not in precalculated_orders or not allow_precalculated:
            # load data from slice
            r_factor_coverage = load_slice_coverage(r_chrm, r_start, r_end)
        else:
            # load precalculated data for entire segment
            r_factor_coverage = load_precalculated_factor_coverage(r_chrm, r_order, r_i)

        # update coverage trackers
        if membership:
            region_membership[:, i] = r_factor_coverage > 0
        else:
            factor_coverage += r_factor_coverage
            region_coverage += (r_end - r_start)
            segment_bin_count += (r_factor_coverage != 0).astype(int)

    if not membership:
        # calculate factor enrichment
        with np.errstate(divide='ignore', invalid='ignore'):
            ooe = (factor_coverage / region_coverage) / gw_factor_coverage
            enr = np.log2(ooe)

    print(f"Time to calculate factor enrichments for region set: {time.time() - start}", flush=True)
    return region_membership if membership else (enr, segment_bin_count) 


"""
Collects the top N factors by enrichment value. Filters out excluded factors. 
"""
def get_top_factors(enr, segment_bin_count, N, factor_exclusion, corank=False, count=True):
    # collect all factors
    factors = load_all_factors()
    
    # exclude factors
    factor_exclusion = np.array([f"{f['dataset']}.{f['factor']}" for f in factor_exclusion])
    exclusion_mask = ~np.isin(factors, factor_exclusion) & (enr > 0)
    factors = factors[exclusion_mask]
    enr = enr[exclusion_mask]
    segment_bin_count = segment_bin_count[exclusion_mask]
    
    if corank:
        enr_corank = rankdata(enr)
        count_corank = rankdata(segment_bin_count)
        top_factor_inds = np.argsort(-(enr_corank + count_corank))[:N]
    elif count:
        # sort by count and take the top N
        top_factor_inds = np.argsort(-segment_bin_count)[:N]
    else:
        # sort by enrichment and take the top N
        top_factor_inds = np.argsort(-enr)[:N]

    return factors[top_factor_inds], enr[top_factor_inds], segment_bin_count[top_factor_inds]


"""
Parses membership for each factor and returns factor names in which region has membership
"""
def parse_membership(region_membership, factor_exclusion=[]):

    start = time.time()

    # load factor names
    factors = load_all_factors()

    # exclude factors
    factor_exclusion = np.array([f"{f['dataset']}.{f['factor']}" for f in factor_exclusion])
    exclusion_mask = ~np.isin(factors, factor_exclusion)
    region_membership = region_membership.flatten()[exclusion_mask]
    factors = factors[exclusion_mask]

    print(f"Time to parse factor membership: {time.time() - start}", flush=True)
    return factors[region_membership]


"""
Tests overlap of factor segments with region. Returns overlapping segments.
"""
def find_overlapping_segments(region, factors_to_test, splitorder=4):

    start = time.time()

    # chromosome of interest
    r_chrm = region['chromosome']
    # splitorder position of region
    r_splitorder_pos = hilbert_pos_to_order(region['i'], region['order'], splitorder)
    h5_dataset_name = f"splitorder_{splitorder}.{r_splitorder_pos}"

    overlapping_segments = {}
    with h5py.File(order_independent_filtering_path, 'r') as f:
        for factor in factors_to_test:

            # load factor data from region split (chromosome and splitorder position)
            h5_group_name = f"{factor}.{r_chrm}"
            factor_data = pickle.loads(f[f"{h5_group_name}/{h5_dataset_name}"][()].tobytes())

            # look at segments at order greater than region order
            above = factor_data[factor_data[:, 0] > region['order']]

            # test overlap
            overlapping = above[
                hilbert_pos_to_order(above[:, 2], above[:, 0], region['order']) == region['i']
            ]

            # track overlapping segments
            if len(overlapping) > 0:
                overlapping_segments[factor] = overlapping

    print(f"Time to find overlapping segments for {len(factors_to_test)} factors: {time.time() - start}", flush=True)
    return overlapping_segments


"""
Refactors factors for return
"""
def generate_overlap_return_data(top_factors, top_enr, top_counts):
    return json.dumps(
        [
            {
                "dataset": f.split('.')[0], "factor": int(f.split('.')[1]), "enrichment": float(e), "count": int(c)
            } for f, e, c in zip(top_factors, top_enr, top_counts)
        ]
    )


"""
Generates return for factor segments overlapping region
"""
def generate_single_region_return_data(overlapping_segments):
    return json.dumps(
        [
            {
                "dataset": f.split('.')[0], "factor": int(f.split('.')[1]), "segments": [
                    {'order': s[0], 'chromosome': s[1], 'i': int(s[2]), 'score': float(s[3])} for s in segments
                ]
            } for f, segments in overlapping_segments.items()
        ]
    )


"""
Endpoint for calculating factor enrichment within a region set
"""
@region_set_enrichment_bp.route("region_set_enrichment", methods=['POST'])
def region_set_enrichment():

    request_info = request.get_json()
    regions = request_info['regions']  # [{chromosome, i, order}, ...]
    N = request_info['N'] if 'N' in request_info else 10  # number of factors to suggest
    factor_exclusion = request_info['factorExclusion'] if 'factorExclusion' in request_info else []  # factors to exclude [{dataset, factor}, ...]
    allow_precalculated = request_info['allow_precalculated'] if 'allow_precalculated' in request_info else True
    allow_precalculated = allow_precalculated.lower() == "true" if type(allow_precalculated) == str else allow_precalculated

    # print("FACTOR EXCLUSION", factor_exclusion, flush=True)

    # test factors for enrichment within region set
    enr, segment_bin_count = get_region_set_factor_enrichments(regions, allow_precalculated=allow_precalculated)
    # collect the top N factors by enrichment
    top_factors, top_enr, top_counts = get_top_factors(enr, segment_bin_count, N, factor_exclusion)
    # refactor for return
    return_arr = generate_overlap_return_data(top_factors, top_enr, top_counts)

    return return_arr


"""
Endpoint for finding single region factor information
"""
@region_set_enrichment_bp.route("single_region_factor_query", methods=['POST'])
def single_region_factor_query():

    request_info = request.get_json()
    region = request_info['region']  # {chromosome, i, order}
    # N = request_info['N'] if 'N' in request_info else 10  # number of factors to suggest
    factor_exclusion = request_info['factorExclusion'] if 'factorExclusion' in request_info else []  # factors to exclude [{dataset, factor}, ...]
    allow_precalculated = request_info['allow_precalculated'] if 'allow_precalculated' in request_info else True
    allow_precalculated = allow_precalculated.lower() == "true" if type(allow_precalculated) == str else allow_precalculated

    # test factors for membership within region
    region_membership = get_region_set_factor_enrichments([region], allow_precalculated=allow_precalculated, membership=True)

    # determine factors for which region has membership
    factors_to_test = parse_membership(region_membership, factor_exclusion)

    # find overlapping segments
    overlapping_segments = find_overlapping_segments(region, factors_to_test)

    return_arr = generate_single_region_return_data(overlapping_segments)

    return return_arr