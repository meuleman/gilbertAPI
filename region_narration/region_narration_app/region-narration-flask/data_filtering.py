from flask import Blueprint, request
import numpy as np
import requests
import json
import h5py
import os
import time
from filtering import fetch_indices_scores

# Create a Blueprint
data_filtering_bp = Blueprint('data_filtering_bp', __name__)

chromosomes = np.append([f"chr{c}" for c in range(1, 23)], ['chrX', 'chrY'])

CHR_SIZES_UCSC = {"chr1": 248956422,
                  "chr10": 133797422,
                  "chr11": 135086622,
                  "chr12": 133275309,
                  "chr13": 114364328,
                  "chr14": 107043718,
                  "chr15": 101991189,
                  "chr16": 90338345,
                  "chr17": 83257441,
                  "chr18": 80373285,
                  "chr19": 58617616,
                  "chr2": 242193529,
                  "chr20": 64444167,
                  "chr21": 46709983,
                  "chr22": 50818468,
                  "chr3": 198295559,
                  "chr4": 190214555,
                  "chr5": 181538259,
                  "chr6": 170805979,
                  "chr7": 159345973,
                  "chr8": 145138636,
                  "chr9": 138394717,
                  "chrX": 156040895,
                  "chrY": 57227415}

datadir="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/data_filtering/20240819"

def hilbert_pos_to_order(pos, frm, to):
    return np.floor(pos / pow(4, frm - to)).astype(int)


"""
For a given chromosome group, find the common indices between all filters.
"""
def find_intersection(group):
    # sort groups by order
    group = [group[i] for i in np.argsort([g['order'] for g in group])]

    indices = group[0]["indices"]
    indices_sort_order = np.argsort(indices)
    indices = indices[indices_sort_order]
    combined_scores = group[0]["scores"][indices_sort_order]

    current_order = group[0]['order']

    # add other groups to the combined indices and scores
    for g in group[1:]:
        g_indices = g["indices"]
        g_sort_order = np.argsort(g_indices)
        g_indices = g_indices[g_sort_order]
        g_scores = g["scores"][g_sort_order]

        # convert current indices at current order to new order
        indices_g_order = hilbert_pos_to_order(indices, current_order, g['order'])

        # find base position for new order indices relative to current order
        g_indices_base = hilbert_pos_to_order(
            hilbert_pos_to_order(
                g_indices, 
                g['order'], 
                current_order
            ), 
            current_order, 
            g['order']
        )

        # find overlap in indices
        common_mask = np.isin(indices_g_order, g_indices_base)
        indices_g_order = indices_g_order[common_mask]
        combined_scores = combined_scores[common_mask]

        # combine scores
        g_mask = np.isin(g_indices_base, indices_g_order)
        g_indices_base = g_indices_base[g_mask]
        g_scores = g_scores[g_mask]

        # count how many times each base position is found with new order indices
        g_indices_base_unique, g_indices_base_counts = np.unique(g_indices_base, return_counts=True)
        # expand the combined scores to match the position occurrences of the new order
        if np.count_nonzero(g_indices_base_unique != indices_g_order) > 0:
            raise ValueError(f"Inconsistent mapping between order {current_order} and order {g['order']} for {g['chromosome']} when calculating intersection")
        combined_scores = np.repeat(combined_scores, g_indices_base_counts)

        indices = g_indices[g_mask]
        combined_scores += g_scores

        # update current order
        current_order = g['order']

    # average combined factor scores
    combined_scores = combined_scores / len(group)
    return {"indices": indices, "scores": combined_scores, 'order': current_order}


def get_filtered_indices_scores(filters):
    filter_indices = []
    start = time.time()
    
    # parse each filter and collect indices
    for f in filters:
        order = f['order']
        index = f['index']
        dataset_name = f['dataset_name']
        # fetch indices from each chromosome
        for chromosome in chromosomes:
            indices_scores = fetch_indices_scores(
                order, 
                chromosome, 
                index, 
                dataset_name, 
                os.path.join(datadir, 'indices.scores.native_resolution.h5')
            )
            if len(indices_scores) > 0:
                filter_indices.append({
                    "chromosome": chromosome, 
                    "order": order, 
                    "index": index, 
                    "dataset_name": dataset_name, 
                    "indices": indices_scores[:, 0].astype(int),
                    "scores": indices_scores[:, 1],
                })

    print(f"Time to fetch indices: {time.time() - start}", flush=True)

    start = time.time()
    # group by chromosome
    chromosome_groups = {}
    for chromosome in chromosomes:
        chromosome_group = [f for f in filter_indices if f['chromosome'] == chromosome]
        # only keep if each filter has indices for this chromosome
        if len(chromosome_group) == len(filters):
            chromosome_groups[chromosome] = chromosome_group


    # reduce chromosome groups to common indices
    reduced_chromosome_groups = {}
    for chromosome, group in chromosome_groups.items():
        group_intersection = find_intersection(group)
        if len(group_intersection["indices"]) > 0:
            reduced_chromosome_groups[chromosome] = group_intersection
    
    print(f"Time to find intersection: {time.time() - start}", flush=True)
    return reduced_chromosome_groups


# """
# Collects the combined factor scores for index positions across chromosomes, sorts by score,
# and returns the indices, chromosomes, and scores for the top N paths across the genome.
# """
# def sort_by_score(chromosome_indices, N):
#     start = time.time()

#     # initialize arrays
#     indices = np.array([], dtype="int32")
#     scores = np.array([], dtype="float32")
#     chromosomes = np.array([])

#     # collect indices and scores from each chromosome
#     for chromosome in chromosome_indices.keys():
#         chrm_indices = chromosome_indices[chromosome]["indices"]
#         chrm_scores = chromosome_indices[chromosome]['scores']
#         indices = np.append(indices, chrm_indices)
#         scores = np.append(scores, chrm_scores)
#         chromosomes = np.append(chromosomes, np.array([chromosome] * len(chrm_indices)))

#     # sort by score
#     sort_order = np.argsort(-scores)[:N]
#     top_scores = scores[sort_order]
#     top_indices = indices[sort_order]
#     top_chromosomes = chromosomes[sort_order]

#     print(f"Time to sort by score: {time.time() - start}", flush=True)

#     return top_indices, top_chromosomes, top_scores


"""
Refactors our top segments into a JSON object for return.
"""
def generate_return_data(chromosome_indices):
    start = time.time()

    # return the filtered indices
    return_indices = []
    for chr in chromosome_indices.keys():
        return_indices.append({
            "chromosome": chr,
            "indices": chromosome_indices[chr]['indices'].tolist(),
            "scores": chromosome_indices[chr]['scores'].tolist(),
            "order": chromosome_indices[chr]['order']
        })

    print(f"Time to generate return data: {time.time() - start}", flush=True)
    return json.dumps({
        "filtered_segments": return_indices
    })


"""
Endpoint for fetching the native order resolution region indices and data for any filter at any order
"""
@data_filtering_bp.route("data_filtering", methods=['POST'])
def data_filtering():
    request_info = request.get_json()
    
    filters = request_info['filters']  # [{order, index, dataset_name}, ...]
    N = request_info['N'] if 'N' in request_info else None  # number of segments to return

    # find the intersection of filters
    chromosome_indices = get_filtered_indices_scores(filters)

    # # find the top segments based on combined factor score
    # top_indices, top_chromosomes, top_scores = sort_by_score(chromosome_indices, N)

    # generate the return list
    return_list = generate_return_data(chromosome_indices)
    
    return return_list