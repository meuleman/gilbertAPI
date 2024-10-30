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
Filters paths down by if there is overlap with any regions. If there is overlap with a 
region and the region is higher resolution than the path, the path is adjusted to include
the region.
"""
def region_filter(paths, regions, min_order):
    start = time.time()

    # group regions by chromosome
    regions_grouped = {}
    for r in regions:
        if r["chromosome"] not in regions_grouped:
            regions_grouped[r["chromosome"]] = []
        regions_grouped[r["chromosome"]].append({"order": r["order"], "i": r["i"]})
    
    confirmed_paths = {}
    for chrm in paths:
        chrm_regions = regions_grouped[chrm] if chrm in regions_grouped else []
        chrm_paths = paths[chrm]

        # seperate orders and positions for paths
        chrm_path_os = np.array([np.where(~np.isnan(p[:, 0]))[0][-1] for p in chrm_paths])
        chrm_path_orders = chrm_path_os + min_order
        chrm_path_positions = np.array([int(p[chrm_path_os[i], 0]) for i, p in enumerate(chrm_paths)])

        # seperate orders and positions for regions
        chrm_region_orders = np.array([r["order"] for r in chrm_regions])
        chrm_region_positions = np.array([r["i"] for r in chrm_regions])

        # for each path, filter by regions
        for i, path_pos, path_order in zip(range(len(chrm_path_positions)), chrm_path_positions, chrm_path_orders):
            r_gt_p_i = np.where(chrm_region_orders > path_order)[0]
            r_lte_p_i = np.where(chrm_region_orders <= path_order)[0]

            # for regions in which order is greater than path order
            # convert regions with order greater than path order to path order
            r_gt_p = chrm_region_positions[r_gt_p_i]
            r_gt_p_orders = chrm_region_orders[r_gt_p_i]
            r_gt_p_converted = hilbert_pos_to_order(r_gt_p, r_gt_p_orders, path_order)

            # get the original regions where the converted regions are matching segment position
            matching_regions = r_gt_p[np.where(r_gt_p_converted == path_pos)[0]]
            matching_region_orders = r_gt_p_orders[np.where(r_gt_p_converted == path_pos)[0]]

            for mr, mro in zip(matching_regions, matching_region_orders):
                if chrm not in confirmed_paths:
                    confirmed_paths[chrm] = []

                # adjust the path to include the region index at the region order
                p = chrm_paths[i].copy()
                p = np.vstack((p, generate_empty_path(length=(mro - min_order + 1 - len(p)))))
                p[mro - min_order] = [mr, 0]  #TODO: incorporate true score?
                confirmed_paths[chrm].append(p)

            # for regions in which order is less than or equal to path order
            r_lte_p = chrm_region_positions[r_lte_p_i]
            r_lte_p_orders = chrm_region_orders[r_lte_p_i]

            # convert path position to region orders
            p_converted = hilbert_pos_to_order(path_pos, path_order, r_lte_p_orders)
            if np.count_nonzero(p_converted == r_lte_p):
                if chrm not in confirmed_paths:
                    confirmed_paths[chrm] = []

                confirmed_paths[chrm].append(chrm_paths[i])

    print(f"Time to filter segments by regions: {time.time() - start}", flush=True)
    return confirmed_paths


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
def get_top_paths(paths, min_order, num_top_paths=None):
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

    while (((num_top_paths == None) or (count < num_top_paths)) and (i < len(winning_scores))):

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
            # test for matches with previous path
            matches = [
                np.count_nonzero(
                    # pad with nans to make same length (in case region was included in path)
                    np.pad(cp[:, 0], (0, max(0, len(p) - len(cp))), constant_values=np.nan) == np.pad(p[:, 0], (0, max(0, len(cp) - len(p))), constant_values=np.nan)
                ) > 0 for cp in top_paths[chrm]
            ]

            # if not positive matches
            if np.count_nonzero(matches) == 0:
                top_paths[chrm].append(p)
                count += 1
            
    print(f"Time to find top paths: {time.time() - start}", flush=True)
    return top_paths


"""
Refactors top paths for return
"""
def generate_return(top_paths, min_order):
    return_array = []
    for chrm in top_paths:

        for p in top_paths[chrm]:
            # max order position in which we see factor
            p_max_o = np.where(~np.isnan(p[:, 0]))[0][-1]

            return_array.append({
                "chromosome": chrm,
                "i": int(p[p_max_o, 0]),
                "order": int(p_max_o + min_order),
                "max_score_order": int(np.nanargmax(p[:, 1]) + min_order),  # not necessarily same as "order"
                "max_score": float(np.nanmax(p[:, 1])),
            })
    print(f"Returning {len(return_array)} regions when filtering without order", flush=True)
    return json.dumps({"regions": return_array})


"""
Endpoint for filtering on a factor without a specified order
"""
@filtering_without_order_bp.route("filtering_without_order", methods=['POST'])
def filtering_without_order():

    request_info = request.get_json()
    factorInd = request_info['factorInd']
    dataset = request_info['dataset']
    regions = request_info['regions'] if 'regions' in request_info else None # [{chromosome, order, i}, ...]

    # construct paths for all segments containing factor
    paths, min_order, max_order = construct_paths(factorInd, dataset)

    # if regions present, make sure paths overlap at least one region
    paths = region_filter(paths, regions) if regions is not None else paths

    # find the top non-overlapping segment paths
    top_paths = get_top_paths(paths, min_order)

    # refactor for return
    return_array = generate_return(top_paths, min_order)

    return return_array