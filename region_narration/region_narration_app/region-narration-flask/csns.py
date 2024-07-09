from flask import Blueprint, request
import numpy as np
import requests
import json
import h5py
import os
import time

# Create a Blueprint
csns_bp = Blueprint('csns_bp', __name__)

baseS3URLPrefix = "https://altius-gilbert.s3.us-west-2.amazonaws.com/20240703/csn_index_files"
chromosomes = np.append([f"chr{c}" for c in range(1, 23)], ['chrX', 'chrY'])

datadir="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/precomputed_csn/20240703"

def hilbert_pos_to_order(pos, frm, to):
    return np.floor(pos / pow(4, frm - to)).astype(int)

"""
Fetches the order 14 positional indices for a given order, chromosome, factor index, 
and dataset name. Returns an empty array if the indices are not found.
"""
def fetch_indices(
    order, 
    chromosome, 
    factor_index, 
    dataset_name, 
    file_path=os.path.join(datadir, 'indices/all_indices.order_14_resolution.h5')
):
    dataset_name = f"{order}.{chromosome}.{dataset_name}.{factor_index}"
    with h5py.File(file_path, 'r') as f:
        try:
            indices = f[dataset_name][:]
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
def get_filtered_indices(filters, filter_chromosomes=[]):
    filter_indices = []
    # parse each filter and collect indices
    start = time.time()
    for f in filters:
        order = f['order']
        index = f['index']
        dataset_name = f['dataset_name']
        # fetch indices from each chromosome
        for chromosome in chromosomes:
            if len(filter_chromosomes) > 0 and chromosome not in filter_chromosomes:
                continue
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
Takes in a dictionary of chromosome indices and filters them by regions.
If more than one regions are present, we OR the regions together, giving all paths that flow through each region
"""
def filter_regions(chromosome_indices, regions):
    filtered_chromosome_indices = {}
    for region in regions:
        chromosome = region['chromosome']
        indices = chromosome_indices[chromosome]
        order = region['order']
        index = region['i']
        ops = hilbert_pos_to_order(indices, 14, order)
        indices = indices[ops == index]
        # we "OR" the regions together if there are more than one (getting all the valid indices in given regions of interest)
        if(filtered_chromosome_indices.get(chromosome)):
            filtered_chromosome_indices[chromosome] = np.append(filtered_chromosome_indices[chromosome], indices)
        else:
            filtered_chromosome_indices[chromosome] = indices
    return filtered_chromosome_indices

"""
Collects the path scores for index positions across chromosome, sorts by score,
and returns the indices, chromosomes, and scores for the top N paths across the genome.
"""
def sort_by_score(chromosome_indices, score_type, N, filter_orders, order_diversity=-1):
    indices = np.array([], dtype="int32")
    chromosomes = np.array([])
    seen_counts = np.array([], dtype="int32")

    start = time.time()
    for chromosome in chromosome_indices.keys():
        chrm_indices = chromosome_indices[chromosome]

        # diversity filter: keep only the first index we find within a segment at highest res filter order
        if order_diversity > 0:
            chrm_ops = hilbert_pos_to_order(chrm_indices, 14, order_diversity)  # same ordering as chrm_indices
            unique_chrm_ops, chrm_ops_indices, chrm_ops_counts = np.unique(
                chrm_ops, return_index=True, return_counts=True
            )  # this sorts the indices

            # we lose the ordering here - we do this to match the count order and get unique original indices
            chrm_indices = chrm_indices[chrm_ops_indices]  
        else:
            chrm_ops_counts = np.ones(len(chrm_indices), dtype="int32")
        
        indices = np.append(indices, chrm_indices)

        chrm_array = np.array([chromosome] * len(chrm_indices))
        chromosomes = np.append(chromosomes, chrm_array)

        seen_counts = np.append(seen_counts, chrm_ops_counts)

    print(f"Time to find diversity: {time.time() - start}", flush=True)

    # get scores
    if ((score_type == "full" )| (score_type == "factor")):
        scores = get_scores(indices, chromosomes, score_type, filter_orders)

    start = time.time()
    # sort
    sort_order = np.argsort(-scores)[:N]
    top_scores = scores[sort_order]
    top_indices = indices[sort_order]
    top_chromosomes = chromosomes[sort_order]
    top_seen_counts = seen_counts[sort_order]

    print(f"Time to sort by score: {time.time() - start}", flush=True)

    return top_indices, top_chromosomes, top_scores, top_seen_counts

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
    score_type,
    filter_orders,
    full_file_path=os.path.join(datadir, 'path_scores/path_scores.order_14_resolution.float32.h5'),
    factor_file_path=os.path.join(datadir, 'factor_scores/factor_scores.order_14_resolution.h5')
):
    start = time.time()
    scores = []
    # load in h5 file and read in the scores for each position
    file_path = full_file_path if score_type == "full" else factor_file_path
    with h5py.File(file_path, 'r') as f:
        unique_chromosomes = np.unique(chromosomes)
        for chromosome in unique_chromosomes:
            chromosome_indices = indices[chromosomes == chromosome]
            sort_order = np.argsort(chromosome_indices)
            unsort_order = np.argsort(sort_order)

            # scores for each segment
            chromosome_scores = f[chromosome][chromosome_indices[sort_order]][unsort_order]
            
            # average scores for factors at each order filtered upon if not path scores
            if score_type == "factor":
                orders_adjusted = np.array(filter_orders) - 4
                chromosome_scores = np.mean(chromosome_scores[:, orders_adjusted], axis=1)

            scores.extend(chromosome_scores)

    print(f"Time to get scores: {time.time() - start}", flush=True)
    return np.array(scores)


"""
Refactors our top paths into a JSON object for return.
"""
def generate_return_data(chromosome_indices, top_indices, top_chromosomes, top_scores, dehydrated_paths, top_seen_counts):
    # return the top csns
    start = time.time()
    return_list = []
    for index, chrom, score, path, seen_count in zip(top_indices, top_chromosomes, top_scores, dehydrated_paths, top_seen_counts):
        return_list.append({
            "i": int(index),
            "chromosome": str(chrom),
            "score": float(score),
            "path": path.tolist(),
            "representedPaths": int(seen_count)
        })

    # return the filtered indices
    # return_indices = []
    # for chr in chromosome_indices.keys():
    #     return_indices.append({
    #         "chromosome": chr,
    #         "indices": chromosome_indices[chr].tolist()
    #     })

    print(f"Time to generate return data: {time.time() - start}", flush=True)
    return json.dumps({
        # "filtered_indices": return_indices,
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

    # figure out which chromosomes the regions are in (if any)
    if len(regions) > 0:
        filter_chromosomes = list(set(r['chromosome'] for r in regions))
    else:
        filter_chromosomes = []
    # get the indices for the filters of each chromosome
    # if we have regions, only get the indices for chromomsomes present in the regions
    chromosome_indices = get_filtered_indices(filters, filter_chromosomes)
    # filter down the indices to the regions, converts indices to order of regions and checks
    chromosome_indices = filter_regions(chromosome_indices, regions)

    filter_orders = [f['order'] for f in filters]
    if(diversity):
        # max_order = max(filters, key=lambda x: x['order'])['order']
        max_order = max(filter_orders)
    else:
        max_order = -1
    # find the top paths based on full or factor score
    # print("sorting by score", flush=True)
    top_indices, top_chromosomes, top_scores, top_seen_counts = sort_by_score(chromosome_indices, score_type, N, filter_orders, order_diversity=max_order)      
    
    # get the dehydrated paths for each index 
    # print("getting dehydrated paths", flush=True)
    dehydrated_paths = get_dehydrated_paths(top_indices, top_chromosomes)

    # generate the return list
    # print("generating return list", flush=True)
    return_list = generate_return_data(chromosome_indices, top_indices, top_chromosomes, top_scores, dehydrated_paths, top_seen_counts)

    return return_list