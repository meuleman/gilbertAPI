from flask import Blueprint, request
from common import hilbert_pos_to_order, fetch_indices_scores
import numpy as np
import h5py
import pickle
# import os
import time
import json
from collections import defaultdict
from constants import chromosomes, PRELOADED_INDICES, orders, filtering_path, order_independent_filtering_path, CHR_SIZES_UCSC

# Create a Blueprint
filtering_without_order_bp = Blueprint('filtering_without_order_bp', __name__)


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
                    chrm_s['scores'].reshape(-1, 1)
                ))
                segments[f"{chrm_s['order']}.{chrm_s['chromosome']}"] = s
                
        # otherwise, fetch from h5 (slower)
        else:
            for chrm in chromosomes:
                s = fetch_indices_scores(order, chrm, factorInd, dataset, filtering_path)
                if len(s) > 0:
                    segments[f"{order}.{chrm}"] = s

    print(f"Time to load filters: {time.time() - start}", flush=True)
    return segments


"""
Collects all segments for a factor and returns them in an array
"""
def get_segments(factor, dataset):
    start = time.time()

    # collect segments
    segments = load_filters(factor, dataset)

    # create segments array
    all_segments = np.array([]).reshape(0, 4)
    for k in list(segments):
        order, chrm = k.split('.')
        order = int(order)

        k_segments = np.hstack((
            np.array([order] * len(segments[k])).reshape(-1, 1),
            np.array([chrm] * len(segments[k])).reshape(-1, 1),
            segments[k]
        ), dtype=object)

        all_segments = np.vstack((all_segments, k_segments))

    # sort segments array first by score, then by order (lowest order first)
    all_segments = all_segments[np.lexsort((all_segments[:, 0], -all_segments[:, 3]))]

    print(f"Time to collect all segments for factor: {time.time() - start}", flush=True)
    return all_segments


"""
A more efficient test_overlap function that uses a bp tracker to keep track of overlapping regions.
Cannot update a segment that overlaps with higher resolution region(s).
"""
def test_overlap(chrm, start, end, tracker, keep_overlapping=False):
    if np.count_nonzero(tracker[chrm][start:end]) == 0:
        return False if keep_overlapping else True
    else:
        return True if keep_overlapping else False
    


"""
Reduces segments to those non-overlapping. Also tests for overlap with background regions if provided.
"""
def test_segments_for_overlap(segments, regions=None):
    start = time.time()

    # find set of non-overlapping segments
    no_segments = []
    overlap_tracker = {
        chrm: np.zeros(CHR_SIZES_UCSC[chrm], dtype=bool) for chrm in np.unique(segments[:, 1])
    }

    # create and fill region tracker
    if regions is not None:
        regions_tracker = {
            chrm: np.zeros(CHR_SIZES_UCSC[chrm], dtype=bool) for chrm in np.unique(segments[:, 1])
        }
        for r in regions:
            r_chrm = r['chromosome']
            if r_chrm in regions_tracker:
                r_start = hilbert_pos_to_order(r['i'], r['order'], 14)
                r_end = hilbert_pos_to_order(r['i'] + 1, r['order'], 14)
                regions_tracker[r_chrm][r_start:r_end] = 1

    for s in segments:
        # bounds of s
        s_chrm = s[1]
        s_start = hilbert_pos_to_order(s[2], s[0], 14)
        s_end = hilbert_pos_to_order(s[2] + 1, s[0], 14)
        add_s = False
        
        # test for overlap
        if test_overlap(s_chrm, s_start, s_end, overlap_tracker):
            if regions is not None:
                if test_overlap(s_chrm, s_start, s_end, regions_tracker, keep_overlapping=True):
                    add_s = True
            else:
                add_s = True

        # if appropriate overlap, add segment and update overlap tracker
        if add_s:
            no_segments.append(s)
            overlap_tracker[s_chrm][s_start:s_end] = 1

    no_segments = np.vstack(no_segments)
    print(f"Time to test segments for overlap: {time.time() - start}", flush=True)
    return no_segments


"""
Refactors segments for return
"""
def generate_return(segments):
    start = time.time()
    return_array = []
    for s in segments:
        return_array.append({
            "chromosome": s[1],
            "i": int(s[2]),
            "order": int(s[0]),
            "score": float(s[3]),
        })

    print(f"Returning {len(return_array)} regions when filtering without order", flush=True)
    print(f"Time to generate return array: {time.time() - start}", flush=True)
    return json.dumps({"regions": return_array})


"""
Looks at segments and background region set and updates segments within an orderXchrm split
to a higher resolution region if overlapping one.
"""
def convert_segment_split(split, against):

    # get the order of split
    order = split[0, 0]

    # find background regions with orders > split order
    above_mask = (against[:, 0] > order)
    above_against = against[above_mask]

    # convert background regions to split order
    above_converted = hilbert_pos_to_order(above_against[:, 2], above_against[:, 0], order)

    # create dictionary to map each region to a split order position
    regions_dict = defaultdict(list)
    for i, pos in enumerate(above_converted):
        regions_dict[pos].append(i)

    # collect regions for each split segment
    segment_updates = [regions_dict[pos] for pos in split[:, 2]]

    # replace segment with region(s) if overlap found
    segments_updated = [
        np.column_stack((
            above_against[u],
            np.full(len(above_against[u]), split[i][3])
        )) if len(u) > 0 else split[i] for i, u in enumerate(segment_updates)
    ]
    return segments_updated


"""
Splits segments into groups by order and chromosome, then within each group updates segments 
to higher resolution regions if any overlap. Returns the sorted updated segments.
"""
def convert_segments(segments, regions):
    start = time.time()

    # collect the unique orders and chromosomes for segments
    segment_orders = np.unique(segments[:, 0])
    segment_chrms = np.unique(segments[:, 1])

    # split segments into groups by order and chromosome
    segments_split = {
        f"{o}.{c}": segments[(segments[:, 0] == o) & (segments[:, 1] == c)] for o in segment_orders for c in segment_chrms
    }
    segments_split = {k: v for k, v in segments_split.items() if len(v) > 0}

    # create region split
    regions_flat = np.array([[r['order'], r['chromosome'], r['i']] for r in regions], dtype=object)
    regions_chrm_split = {
        chrm: regions_flat[regions_flat[:, 1] == chrm] if chrm in segment_chrms else None for chrm in chromosomes
    }

    # for each split, match background regions and see if segments need to be updated
    segments_updated = []
    for split in list(segments_split.values()):
        chrm = split[0, 1]
        segments_updated.extend(convert_segments(split, regions_chrm_split[chrm]))

    segments_updated = np.vstack(segments_updated)
    # sort first by score, then by order (lowest order first)
    segments_updated = segments_updated[np.lexsort((segments_updated[:, 0], -segments_updated[:, 3]))]

    print(f"Time to update segments to higher resolution: {time.time() - start}", flush=True)
    return segments_updated


"""
Loads precalculated, non-overlapping segments for factor and dataset
"""
def load_precalculated_segments(order_independent_filtering_path, factor, dataset):
    start = time.time()

    # load precalculated segments for factor and dataset
    with h5py.File(order_independent_filtering_path, 'r') as f:
        segments = pickle.loads(f[f"{dataset}.{factor}"][()].tobytes())
        
    print(f"Time to load precalculated segments: {time.time() - start}")
    return segments


"""
Endpoint for filtering on a factor without a specified order
"""
@filtering_without_order_bp.route("filtering_without_order", methods=['POST'])
def filtering_without_order():

    request_info = request.get_json()
    factorInd = request_info['factorInd']
    dataset = request_info['dataset']
    regions = request_info['regions'] if 'regions' in request_info else None # [{chromosome, order, i}, ...]
    force_recalculation = request_info['force_recalculation'] if 'force_recalculation' in request_info else False
    force_recalculation = force_recalculation.lower() == "true" if type(force_recalculation) == str else force_recalculation

    # get segments and test for overlap
    if force_recalculation or regions is not None:
        # get segments
        segments = get_segments(factorInd, dataset)

        # test segments for overlap with each other and with background regions
        no_segments = test_segments_for_overlap(segments, regions)

        # convert segments to higher resolution regions
        if regions is not None:
            no_segments = convert_segments(no_segments, regions)
    
    # load precalculated, non-overlapping segments
    else:
        no_segments = load_precalculated_segments(order_independent_filtering_path, factorInd, dataset)

    # refactor for return
    return_array = generate_return(no_segments)

    return return_array