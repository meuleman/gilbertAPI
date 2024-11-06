from flask import Blueprint, request
from common import hilbert_pos_to_order, fetch_indices_scores
import numpy as np
# import h5py
# import os
import time
import json
from constants import chromosomes, PRELOADED_INDICES, orders, filtering_path, CHR_SIZES_UCSC

# Create a Blueprint
filtering_without_order_bp = Blueprint('filtering_without_order_bp', __name__)


# """
# Function to load path for a position
# """
# def load_path(chrm, position, dataset='factor'):
    
#     path = os.path.join(datadir, 'paths/full_path_factors.order_14_resolution.h5') if dataset == 'factor' else os.path.join(datadir, 'scores/factor_scores.order_14_resolution.h5')
    
#     with h5py.File(path, 'r') as f:
#         path_data = f[chrm][position]
    
#     return path_data


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

    # sort segments array
    all_segments = all_segments[np.argsort(-all_segments[:, 3])]

    print(f"Time to collect all segments for factor: {time.time() - start}", flush=True)
    return all_segments


# """
# For a given segment, tests overlap with background region set. Region set must be in the format:
# np.array([[order, chromosome, position], ...])
# If keep_overlapping is True, returns True if there is overlap. If False, returns True if there is no overlap.
# If return_s is True, returns the segment and if the segment overlaps a higher resolution region (or multiple higher resolution regions), 
# replaces the segment with the higher resolution region(s).
# """
# def test_overlap(s, against, keep_overlapping=False, return_s=False):
#     below_mask = (against[:, 0] < s[0]) & (against[:, 1] == s[1])
#     above_mask = (against[:, 0] >= s[0]) & (against[:, 1] == s[1])

#     # below
#     s_converted = hilbert_pos_to_order(s[2], s[0], against[below_mask][:, 0])
#     below_count = np.count_nonzero(against[below_mask][:, 2] == s_converted)

#     # above
#     above_converted = hilbert_pos_to_order(against[above_mask][:, 2], against[above_mask][:, 0], s[0])
#     above_count = np.count_nonzero(above_converted == s[2])
    
#     overlap = ((below_count > 0) or (above_count > 0)) if keep_overlapping else ((below_count == 0) and (above_count == 0))
    
#     # replace s with higher resolution region
#     if return_s:
#         if above_count > 0:
#             above_matching = against[above_mask][above_converted == s[2]]
#             s = np.hstack((
#                 above_matching,
#                 np.array([s[3]] * len(above_matching)).reshape(-1, 1)
#             ))
#         else:
#             s = np.array([s])

#     return (overlap, s) if return_s else overlap


# """
# Reduces segments to those non-overlapping. Also tests for overlap with background regions if provided.
# """
# def test_segments_for_overlap(segments, regions=None):
#     start = time.time()

#     # flatten regions
#     regions_flat = np.array([[r['order'], r['chromosome'], r['i']] for r in regions], dtype=object) if regions is not None else []
    
#     # initialize non overlapping array
#     no_segments = np.array([], dtype=object).reshape(0, segments.shape[1])

#     for s in segments:
#         # make sure segment does not overlap previously added segments
#         if test_overlap(s, no_segments):
#             # if background regions provided
#             if len(regions_flat) > 0:
#                 # make sure segment overlaps at least one region
#                 region_overlap, s_arr = test_overlap(s, regions_flat, keep_overlapping=True, return_s=True)
#                 if region_overlap:
#                     no_segments = np.append(no_segments, s_arr, axis=0)
#             else:
#                 no_segments = np.append(no_segments, s.reshape(1, -1), axis=0)

#     print(f"Time to test segments for overlap: {time.time() - start}", flush=True)
#     return no_segments


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
    no_segments = np.array([], dtype=object).reshape(0, segments.shape[1])
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
            r_start = hilbert_pos_to_order(r['i'], r['order'], 14)
            r_end = hilbert_pos_to_order(r['i'] + 1, r['order'], 14)
            regions_tracker[r_chrm][r_start:r_end] = 1

    for s in segments:
        # bounds of s
        s_chrm = s[1]
        s_start = hilbert_pos_to_order(s[2], s[0], 14)
        s_end = hilbert_pos_to_order(s[2] + 1, s[0], 14)
        
        # test for overlap
        if test_overlap(s_chrm, s_start, s_end, overlap_tracker):
            if regions is not None:
                if test_overlap(s_chrm, s_start, s_end, regions_tracker, keep_overlapping=True):
                    no_segments = np.append(no_segments, s.reshape(1, -1), axis=0)
                    overlap_tracker[s_chrm][s_start:s_end] = 1
                    
            else:
                no_segments = np.append(no_segments, s.reshape(1, -1), axis=0)
                overlap_tracker[s_chrm][s_start:s_end] = 1

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
Endpoint for filtering on a factor without a specified order
"""
@filtering_without_order_bp.route("filtering_without_order", methods=['POST'])
def filtering_without_order():

    request_info = request.get_json()
    factorInd = request_info['factorInd']
    dataset = request_info['dataset']
    regions = request_info['regions'] if 'regions' in request_info else None # [{chromosome, order, i}, ...]

    # get segments
    segments = get_segments(factorInd, dataset)

    # test segments for overlap with each other and with background regions
    no_segments = test_segments_for_overlap(segments, regions)

    # refactor for return
    return_array = generate_return(no_segments)

    return return_array