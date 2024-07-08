from flask import Blueprint, request
import numpy as np
import requests
import json
import h5py
import os
import time
from math import floor

# Create a Blueprint
csns_bp = Blueprint('csns_bp', __name__)

baseS3URLPrefix = "https://altius-gilbert.s3.us-west-2.amazonaws.com/20240703/csn_index_files"
chromosomes = np.append([f"chr{c}" for c in range(1, 23)], ['chrX', 'chrY'])

datadir="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/precomputed_csn/20240703"

def hilbert_pos_to_order(pos, frm, to):
    return floor(pos / pow(4, frm - to))

"""
Fetches the order 14 positional indices for a given order, chromosome, factor index, 
and dataset name. Returns an empty array if the indices are not found.
"""
def fetch_indices(order, chromosome, factor_index, dataset_name):
    url = f"{baseS3URLPrefix}/{order}.{chromosome}.{dataset_name}.{factor_index}.order_14_resolution.indices.int32.bytes"
    try:
        bytes_arr = requests.get(url).content
        indices = np.frombuffer(bytes_arr, dtype="int32")
    except:
        indices = np.array([], dtype="int32")
    return indices


"""
For a given chromosome group, find the common indices between all filters.
"""
def find_intersection(group):
    common_indices = group[0]["indices"]
    for g in group[1:]:
        common_indices = common_indices[np.isin(common_indices, g["indices"])]
    return common_indices


"""
Collects positional indices for each filter, groups them by chromosome, and finds the
intersection between filters within each chromosome.
"""
def get_filtered_indices(filters):
    filter_indices = []
    # parse each filter and collect indices
    start = time.time()
    for f in filters:
        order = f['order']
        index = f['index']
        dataset_name = f['dataset_name']
        # fetch indices from each chromosome
        for chromosome in chromosomes:
            indices = fetch_indices(order, chromosome, index, dataset_name)
            if len(indices) > 0:
                filter_indices.append({
                    "chromosome": chromosome, 
                    "order": order, 
                    "index": index, 
                    "dataset_name": dataset_name, 
                    "indices": indices
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
        if len(group_intersection) > 0:
            reduced_chromosome_groups[chromosome] = group_intersection
    
    print(f"Time to find intersection: {time.time() - start}", flush=True)
    return reduced_chromosome_groups


"""
Collects the path scores for index positions across chromosome, sorts by score,
and returns the indices, chromosomes, and scores for the top N paths across the genome.
"""
def sort_by_score(chromosome_indices, score_type, N, order_diversity = -1):
    indices = np.array([], dtype="int32")
    chromosomes = []

    for chromosome in chromosome_indices.keys():
        chrm_indices = chromosome_indices[chromosome]
        chrm_array = np.array([chromosome] * len(chrm_indices))
        
        indices = np.append(indices, chrm_indices)
        chromosomes = np.append(chromosomes, chrm_array)

    # diversity filter: keep only the first index we find within a segment at highest res filter order
    if order_diversity > 0:
        start = time.time()
        # we use the highest resolution order in the filter
        seen = {}
        diverse_indices = []
        diverse_chromosomes = []
        for index in indices:
            order = hilbert_pos_to_order(index, 14, order_diversity)
            if order not in seen:
                seen[order] = index
                diverse_indices.append(index)
                diverse_chromosomes.append(chromosomes[index])

        print("len all indices", len(indices))
        print("len diverse indices", len(diverse_indices))
        indices = np.array(diverse_indices, dtype="int32")
        chromosomes = np.array(diverse_chromosomes, dtype="str")
        print(f"Time to find diversity: {time.time() - start}", flush=True)

    # get scores
    if score_type == "full":
        scores = get_scores(indices, chromosomes)

    start = time.time()
    # sort
    sort_order = np.argsort(-scores)[:N]
    top_scores = scores[sort_order]
    top_indices = indices[sort_order]
    top_chromosomes = chromosomes[sort_order]

    print(f"Time to sort by score: {time.time() - start}", flush=True)

    return top_indices, top_chromosomes, top_scores

"""
Collects and returns the factor indices for each of the top paths.
"""
def get_dehydrated_paths(
    top_indices,
    top_chromosomes,
    file_path=os.path.join(datadir, 'paths/path_factor_indices.order_14_resolution.h5')
):
    start = time.time()
    dehydrated_paths = []
    # load in h5 file and read in the exact dehydrated paths we need
    with h5py.File(file_path, 'r') as f:
        for chromosome, index in zip(top_chromosomes, top_indices):
            dehydrated_paths.append(f[chromosome][index])
    
    print(f"Time to get dehydrated paths: {time.time() - start}", flush=True)
    return dehydrated_paths

"""
Collects path scores for index positions on each chromosome.
"""
def get_scores(
    indices,
    chromosomes,
    file_path=os.path.join(datadir, 'path_scores/path_scores.order_14_resolution.float32.h5')
):
    start = time.time()
    scores = []
    # load in h5 file and read in the scores for each position
    with h5py.File(file_path, 'r') as f:
        unique_chromosomes = np.unique(chromosomes)
        for chromosome in unique_chromosomes:
            chromosome_indices = indices[chromosomes == chromosome]
            sort_order = np.argsort(chromosome_indices)
            unsort_order = np.argsort(sort_order)
            chromosome_scores = f[chromosome][chromosome_indices[sort_order]][unsort_order]
            scores.extend(chromosome_scores)

    print(f"Time to get scores: {time.time() - start}", flush=True)
    return np.array(scores)


"""
Refactors our top paths into a JSON object for return.
"""
def generate_return_data(chromosome_indices, top_indices, top_chromosomes, top_scores, dehydrated_paths):
    # return the top csns
    start = time.time()
    return_list = []
    for index, chrom, score, path in zip(top_indices, top_chromosomes, top_scores, dehydrated_paths):
        return_list.append({
            "i": int(index),
            "chromosome": str(chrom),
            "score": float(score),
            "path": path.tolist()
        })

    # return the filtered indices
    return_indices = []
    for chr in chromosome_indices.keys():
        return_indices.append({
            "chromosome": chr,
            "indices": chromosome_indices[chr].tolist()
        })

    print(f"Time to generate return data: {time.time() - start}", flush=True)
    return json.dumps({
        "filtered_indices": return_indices,
        "csns": return_list
    })


"""
Endpoint for fetching the top paths based on a set of filters and regions.
"""
@csns_bp.route("top_paths", methods=['POST'])
def top_paths():
    request_info = request.get_json()
    
    filters = request_info['filters']  # [{order, index}, ...]
    regions = request_info['regions']  # [{order, chromosome, index},...]
    score_type = request_info['scoreType']  # "full", "factor"
    diversity = request_info['diversity']  # true | false
    N = request_info['N']  # number of paths to return

    # get the indices for the filters of each chromosome
    # print("getting filtered indices", flush=True)
    chromosome_indices = get_filtered_indices(filters)

    #TODO: filter by regions that are passed in


    if(diversity):
        max_order = max(filters, key=lambda x: x['order'])['order']
    else:
        max_order = -1
    # find the top paths based on full or factor score
    # print("sorting by score", flush=True)
    top_indices, top_chromosomes, top_scores = sort_by_score(chromosome_indices, score_type, N, order_diversity=max_order)      
    
    # get the dehydrated paths for each index 
    # print("getting dehydrated paths", flush=True)
    dehydrated_paths = get_dehydrated_paths(top_indices, top_chromosomes)

    # generate the return list
    # print("generating return list", flush=True)
    return_list = generate_return_data(chromosome_indices, top_indices, top_chromosomes, top_scores, dehydrated_paths)

    return return_list