from flask import Blueprint, request
from common import hilbert_pos_to_order
import numpy as np
import h5py
# import pickle
# import os
import time
import json
from constants import order_independent_filtering_path

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
Calculates the enrichment of bp coverage for each factor within a region set. Returns 
the log2 of the observed over expected coverage for each factor.
"""
def get_region_set_factor_enrichments(regions):
    start = time.time()
    # load the genome-wide coverage of each factor
    gw_factor_coverage = load_gw_coverage()

    # for each region in region set, determine coverage of each factor
    region_coverage = 0
    factor_coverage = np.zeros(len(gw_factor_coverage))
    for r in regions:
        r_chrm = r['chromosome']
        r_start = hilbert_pos_to_order(r["i"], r["order"], 14)
        r_end = hilbert_pos_to_order(r["i"] + 1, r["order"], 14)

        # collect coverage for each factor within region
        r_factor_coverage = load_slice_coverage(r_chrm, r_start, r_end)

        # update coverage trackers
        factor_coverage += r_factor_coverage
        region_coverage += (r_end - r_start)

    # calculate factor enrichment
    with np.errstate(divide='ignore', invalid='ignore'):
        ooe = (factor_coverage / region_coverage) / gw_factor_coverage
        enr = np.log2(ooe)

    print(f"Time to calculate factor enrichments for region set: {time.time() - start}", flush=True)
    return enr


"""
Collects the top N factors by enrichment value. Filters out excluded factors. 
"""
def get_top_factors(enr, N, factor_exclusion):
    # collect all factors
    factors = load_all_factors()
    
    # exclude factors
    factor_exclusion = np.array([f"{f['dataset']}.{f['index']}" for f in factor_exclusion])
    exclusion_mask = ~np.isin(factors, factor_exclusion)
    factors = factors[exclusion_mask]
    enr = enr[exclusion_mask]
    
    # sort by enrichment
    sort_order = np.argsort(-enr)
    # choose top N (ensuring positive enrichment)
    top_factor_inds = sort_order[enr[sort_order] > 0][:N]

    return factors[top_factor_inds], enr[top_factor_inds]


"""
Refactors factors for return
"""
def generate_return_data(top_factors, top_enr):
    return json.dumps([{"dataset": f.split('.')[0], "factor": int(f.split('.')[1]), "enrichment": e} for f, e in zip(top_factors, top_enr)])


"""
Endpoint for calculating factor enrichment within a region set
"""
@region_set_enrichment_bp.route("region_set_enrichment", methods=['POST'])
def filtering_without_order():

    request_info = request.get_json()
    regions = request_info['regions']  # [{chromosome, i, order}, ...]
    N = request_info['N'] if 'N' in request_info else 10  # number of factors to suggest
    factor_exclusion = request_info['factor_exclusion'] if 'factor_exclusion' in request_info else []  # factors to exclude [{dataset, index}, ...]

    # test factors for enrichment within region set
    enr = get_region_set_factor_enrichments(regions)
    # collect the top N factors by enrichment
    top_factors, top_enr = get_top_factors(enr, N, factor_exclusion)
    # refactor for return
    return_arr = generate_return_data(top_factors, top_enr)

    return return_arr