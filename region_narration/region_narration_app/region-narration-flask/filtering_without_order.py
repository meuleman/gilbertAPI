from flask import Blueprint, request
from common import hilbert_pos_to_order, fetch_indices_scores
import numpy as np
import h5py
import os
import time
import json
from constants import chromosomes, PRELOADED_INDICES, orders, csn_dir as datadir, filtering_path

# Create a Blueprint
filtering_without_order_bp = Blueprint('filtering_without_order_bp', __name__)


"""
Function to load path for a position
"""
def load_path(chrm, position, dataset='factor'):
    
    path = os.path.join(datadir, 'paths/full_path_factors.order_14_resolution.h5') if dataset == 'factor' else os.path.join(datadir, 'scores/factor_scores.order_14_resolution.h5')
    
    with h5py.File(path, 'r') as f:
        path_data = f[chrm][position]
    
    return path_data


"""
Function to load all filters for a single factor
"""
def load_filters(factorInd, dataset, force_h5=False):
    start = time.time()
    
    segments = {}
    for order in orders:
        # if not force_h5, fetch preloaded indices
        if ((not force_h5) & (int(order) <= PRELOADED_INDICES['max_order'])):
            key = f"{order}.{dataset}.{factorInd}"
            s_across_chrm = PRELOADED_INDICES[key] if key in PRELOADED_INDICES else []
            for chrm_s in s_across_chrm:
                s = np.hstack((
                    chrm_s['indices'].reshape(-1, 1),
                    chrm_s['scores'].reshape(-1, 1),
                    np.zeros(len(chrm_s['indices'])).reshape(-1, 1)
                ))
                segments[f"{chrm_s['order']}.{chrm_s['chromosome']}"] = s
                
        # otherwise, fetch from h5 (slower)
        else:
            for chrm in chromosomes:
                s = fetch_indices_scores(order, chrm, factorInd, dataset, filtering_path)
                if len(s) > 0:
                    s = np.hstack((s, np.zeros(len(s)).reshape(-1, 1)))
                    segments[f"{order}.{chrm}"] = s

    print(f"Time to load filters: {time.time() - start}", flush=True)
    return segments


"""
Generates an empty path of specified length
"""
def generate_empty_path(length=11):
    return np.zeros((length, 2)) * np.nan


"""
Collects all segments containing the factor and creates paths from them
"""
def construct_paths(factorInd, dataset):
    start = time.time()

    # collect all filters for factor
    segments = load_filters(factorInd, dataset)
    segments_dropdown = segments.copy()

    # what is the max order with results
    segment_keys = list(segments.keys())
    max_order = max([int(k.split('.')[0]) for k in segment_keys])
    min_order = 4
    max_o = max_order - min_order

    # start search with the highest orders first as lower order segments 
    # can occur within multiple higher order segments
    paths = {}
    for k in np.flip(segment_keys):
        order, chrm = k.split('.')
        order = int(order)
        o = order - min_order

        # add chrm if not yet represented
        if chrm not in paths:
            paths[chrm] = []

        chrm_paths = paths[chrm]

        # first go through current paths and add what is available
        for p in chrm_paths:
            p_o = np.where(~np.isnan(p[:, 0]))[0][-1]
            p_i = p[p_o, 0]
            p_i_converted = hilbert_pos_to_order(p_i, min_order + p_o, order)
            matching = np.where(segments_dropdown[k][:, 0] == p_i_converted)[0]
            if len(matching) == 1:
                p[o] = segments_dropdown[k][matching[0]][:2]
                segments_dropdown[k][matching[0]][-1] = 1  # mark as added

        # for all unadded segments at order, start a new path and add
        for s in segments_dropdown[k][segments_dropdown[k][:, 2] != 1]:
            new_path = generate_empty_path(length=max_o + 1)
            new_path[o] = s[:2]
            chrm_paths.append(new_path)

    print(f"Time to construct segment paths: {time.time() - start}", flush=True)
    
    return paths, min_order, max_order


"""
Function to find the top non-overlapping paths for a factor.
"""
def get_top_paths(paths, min_order, num_top_paths=100):
    start = time.time()

    # chromosomes with at least 1 path
    path_chrms = np.array(list(paths.keys()))
    # max factor scores in each path
    winning_scores = np.array([np.nanmax(p[:, 1]) for k in path_chrms for p in paths[k]])
    # sort scores
    top_winning_score_inds = np.argsort(-winning_scores)
    # the number of paths for each chromosome
    chrm_lens = np.array([len(paths[k]) for k in path_chrms])
    # in order, cumulative sum of all results across chromosomes
    chrm_lens_cumsum = np.cumsum(chrm_lens)

    # add paths if not overlapping previously added
    top_paths = {}
    count = 0
    i = 0

    while ((count < num_top_paths) & (i < len(winning_scores))):

        ind = top_winning_score_inds[i]  # sorted by score
        i += 1

        # determine chromosome and position in chromosome results
        c = np.count_nonzero(chrm_lens_cumsum <= ind)
        chrm = path_chrms[c]
        chrm_ind = ind - (chrm_lens_cumsum[c - 1] if c > 0 else 0)

        # path
        p = paths[chrm][chrm_ind]

        # test for overlap with previously added paths
        if chrm not in top_paths:
            top_paths[chrm] = [p]
            count += 1
        else:
            chrm_paths = top_paths[chrm]
            # test for matches with previous path
            matches = [np.count_nonzero(cp[:, 0] == p[:, 0]) > 0 for cp in chrm_paths]

            # if not positive matches
            if np.count_nonzero(matches) == 0:
                top_paths[chrm].append(p)
                count += 1
            
    print(f"Time to find top paths: {time.time() - start}", flush=True)
    return top_paths


"""
Collect CSN paths for each provided top_path constructed from segments
"""
def get_csn_paths(top_paths, min_order, to_order=14):
    start = time.time()

    top_csn_paths = {}
    top_csn_path_scores = {}
    # each key is a chromosome
    for k in list(top_paths.keys()):

        top_csn_paths[k] = []
        top_csn_path_scores[k] = []
        # each path within a chromosome
        for p in top_paths[k]:
            # max order position in which we see factor
            p_max_o = np.where(~np.isnan(p[:, 0]))[0][-1]
            # convert to order
            p_max_order = p_max_o + min_order
            # position at highest order
            p_order_ind = int(p[p_max_o, 0])
            # convert to order 14
            p_14_ind = hilbert_pos_to_order(p_order_ind, p_max_order, to_order)

            # get the true path and scores
            true_p = load_path(k, p_14_ind, dataset='factor')
            true_p_scores = load_path(k, p_14_ind, dataset='score')

            top_csn_paths[k].append(true_p)
            top_csn_path_scores[k].append(true_p_scores)
            
    print(f"Time to collect CSN paths for segment paths: {time.time() - start}", flush=True)
    return top_csn_paths, top_csn_path_scores


"""
Refactors top paths for return
"""
def generate_return(top_paths, top_csn_paths, top_csn_path_scores, min_order):
    return_array = []
    for chrm in top_paths:
        chrm_paths = top_paths[chrm]
        chrm_csn_paths = top_csn_paths[chrm]
        chrm_csn_path_scores = top_csn_path_scores[chrm]

        for p, p_csn, p_scores in zip(chrm_paths, chrm_csn_paths, chrm_csn_path_scores):
            # max order position in which we see factor
            p_max_o = np.where(~np.isnan(p[:, 0]))[0][-1]

            return_array.append({
                "chromosome": chrm,
                "i": int(p[p_max_o, 0]),
                "order": int(p_max_o + min_order),
                "max_score_order": int(np.nanargmax(p[:, 1]) + min_order),
                "csn": p_csn.tolist(),
                "csn_scores": p_scores.tolist()
            })

    return json.dumps({"top_paths": return_array})


"""
Endpoint for filtering on a factor without a specified order
"""
@filtering_without_order_bp.route("filtering_without_order", methods=['POST'])
def filtering_without_order():

    request_info = request.get_json()
    factorInd = request_info['factorInd']
    dataset = request_info['dataset']

    # construct paths for all segments containing factor
    paths, min_order, max_order = construct_paths(factorInd, dataset)
    # find the top non-overlapping segment paths
    top_paths = get_top_paths(paths, min_order)
    # collect CSN paths for each top segment path
    top_csn_paths, top_csn_path_scores = get_csn_paths(top_paths, min_order)
    # refactor for return
    return_array = generate_return(top_paths, top_csn_paths, top_csn_path_scores, min_order)

    return return_array