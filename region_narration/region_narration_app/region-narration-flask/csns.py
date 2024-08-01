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

datadir="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/precomputed_csn/20240717"

def hilbert_pos_to_order(pos, frm, to):
    return np.floor(pos / pow(4, frm - to)).astype(int)

"""
Fetches the order 14 positional indices for a given order, chromosome, factor index, 
and dataset name. Returns an empty array if the indices are not found.
"""
def fetch_indices_scores(
    order, 
    chromosome, 
    factor_index, 
    dataset_name, 
    file_path=os.path.join(datadir, 'indices/indices.scores.order_14_resolution.h5')
):
    dataset_name = f"{order}.{chromosome}.{dataset_name}.{factor_index}"
    with h5py.File(file_path, 'r') as f:
        try:
            indices_scores = f[dataset_name][:]
        except:
            indices_scores = np.array([], dtype="int32")
    return indices_scores


"""
For a given chromosome group, find the common indices between all filters.
"""
def find_intersection(group):
    indices = group[0]["indices"]
    indices_sort_order = np.argsort(indices)
    indices = indices[indices_sort_order]
    combined_factor_scores = group[0]["factor_scores"][indices_sort_order]
    path_scores = group[0]["path_scores"][indices_sort_order]
    
    # add other groups to the combined indices and scores
    for g in group[1:]:
        g_indices = g["indices"]
        g_sort_order = np.argsort(g_indices)
        g_indices = g_indices[g_sort_order]
        g_factor_scores = g["factor_scores"][g_sort_order]
        
        # find overlap in indices
        common_mask = np.isin(indices, g_indices)
        indices = indices[common_mask]
        path_scores = path_scores[common_mask]
        combined_factor_scores = combined_factor_scores[common_mask]
       
       # combine factor scores
        g_mask = np.isin(g_indices, indices)
        combined_factor_scores += g_factor_scores[g_mask]
    
    # average combined factor scores
    combined_factor_scores = combined_factor_scores / len(group)
    return {"indices": indices, "factor_scores": combined_factor_scores, "path_scores": path_scores}


"""
Collects positional indices for each filter, groups them by chromosome, and finds the
intersection between filters within each chromosome.
"""
def get_filtered_indices_and_scores(filters, region):
    filter_indices = []
    start = time.time()

    # if we have no filters but we do have a region
    # we want to use the range(start, end) in order 14 of the given region
    if len(filters) == 0 and region is not None:
        print("filtering a region with no filters", region)
        chromosome = region["chromosome"]
        start = hilbert_pos_to_order(region["i"], region["order"], 14)
        end = hilbert_pos_to_order(region["i"] + 1, region["order"], 14)
        end = min(end, CHR_SIZES_UCSC[chromosome])
        indices = np.arange(start, end)
        print("len indices", len(indices))
        # no need to deal with other chromosomes
        return {
            chromosome: indices
        }

    # parse each filter and collect indices
    for f in filters:
        order = f['order']
        index = f['index']
        dataset_name = f['dataset_name']
        # fetch indices from each chromosome
        for chromosome in chromosomes:
            if region is not None and chromosome != region["chromosome"]:
                continue
            indices_scores = fetch_indices_scores(order, chromosome, index, dataset_name)
            if len(indices_scores) > 0:
                filter_indices.append({
                    "chromosome": chromosome, 
                    "order": order, 
                    "index": index, 
                    "dataset_name": dataset_name, 
                    "indices": indices_scores[:, 0].astype(int),
                    "path_scores": indices_scores[:, 1],
                    "factor_scores": indices_scores[:, 2]
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
    # print(("CHROMOSOME GROUPS:", chromosome_groups), flush=True)
    for chromosome, group in chromosome_groups.items():
        print("GROUP", group, flush=True)
        group_intersection = find_intersection(group)
        print("GROUP INTERSECTION", group_intersection, flush=True)
        if len(group_intersection["indices"]) > 0:
            reduced_chromosome_groups[chromosome] = group_intersection
    
    print(f"Time to find intersection: {time.time() - start}", flush=True)
    return reduced_chromosome_groups

"""
Takes in a dictionary of chromosome indices and filters them by regions.
If more than one regions are present, we OR the regions together, giving all paths that flow through each region
"""
def filter_region(chromosome_indices, region):
    filtered_chromosome_indices = {}
    # for region in regions:
    chromosome = region['chromosome']
    if chromosome not in chromosome_indices:
        filtered_chromosome_indices[chromosome] = {
            "indices": np.array([], dtype="int32"),
            "path_scores": np.array([], dtype="float32"),
            "factor_scores": np.array([], dtype="float32")
        }
        return filtered_chromosome_indices
    indices = chromosome_indices[chromosome]["indices"]
    path_scores = chromosome_indices[chromosome]["path_scores"]
    factor_scores = chromosome_indices[chromosome]["factor_scores"]
    order = region['order']
    index = region['i']
    ops = hilbert_pos_to_order(indices, 14, order)
    mask = ops == index
    indices = indices[mask]
    path_scores = path_scores[mask]
    factor_scores = factor_scores[mask]
    # we "OR" the regions together if there are more than one (getting all the valid indices in given regions of interest)
    # if(filtered_chromosome_indices.get(chromosome)):
    #     filtered_chromosome_indices[chromosome] = np.append(filtered_chromosome_indices[chromosome], indices)
    # else:
    filtered_chromosome_indices[chromosome] = {"indices": indices, "path_scores": path_scores, "factor_scores": factor_scores}
    return filtered_chromosome_indices

"""
Collects the path scores for index positions across chromosome, sorts by score,
and returns the indices, chromosomes, and scores for the top N paths across the genome.
"""
def sort_by_score(chromosome_indices, score_type, N, filter_orders, order_diversity=-1):
    indices = np.array([], dtype="int32")
    path_scores = np.array([], dtype="float32")
    factor_scores = np.array([], dtype="float32")
    chromosomes = np.array([])
    seen_counts = np.array([], dtype="int32")

    start = time.time()
    for chromosome in chromosome_indices.keys():
        chrm_indices = chromosome_indices[chromosome]["indices"]
        chrm_factor_scores = chromosome_indices[chromosome]["factor_scores"]
        chrm_path_scores = chromosome_indices[chromosome]["path_scores"]

        # diversity filter: keep only the first index we find within a segment at highest res filter order
        if order_diversity > 0:
            chrm_ops = hilbert_pos_to_order(chrm_indices, 14, order_diversity)  # same ordering as chrm_indices
            unique_chrm_ops, chrm_ops_indices, chrm_ops_counts = np.unique(
                chrm_ops, return_index=True, return_counts=True
            )  # this sorts the indices

            # we lose the ordering here - we do this to match the count order and get unique original indices
            chrm_indices = chrm_indices[chrm_ops_indices] 
            chrm_factor_scores = chrm_factor_scores[chrm_ops_indices]
            chrm_path_scores = chrm_path_scores[chrm_ops_indices]
        else:
            chrm_ops_counts = np.ones(len(chrm_indices), dtype="int32")
        
        indices = np.append(indices, chrm_indices)
        path_scores = np.append(path_scores, chrm_path_scores)
        factor_scores = np.append(factor_scores, chrm_factor_scores)
        chromosomes = np.append(chromosomes, np.array([chromosome] * len(chrm_indices)))
        seen_counts = np.append(seen_counts, chrm_ops_counts)

    print(f"Time to find diversity: {time.time() - start}", flush=True)

    # get scores
    # if ((score_type == "full" )| (score_type == "factor")):
    #     scores = get_scores(indices, chromosomes, score_type, filter_orders)
    if (score_type == "full"):
        scores = path_scores
    elif (score_type == "factor"):
        scores = factor_scores

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
    file_path=os.path.join(datadir, 'paths/full_path_factors.order_14_resolution.h5')
):
    start = time.time()
    dehydrated_paths = []
    # load in h5 file and read in the exact dehydrated paths we need
    with h5py.File(file_path, 'r') as f:
        for chromosome, index in zip(top_chromosomes, top_indices):
            dehydrated_paths.append(f[chromosome][index])
    
    print(f"Time to get dehydrated paths: {time.time() - start}", flush=True)
    return dehydrated_paths

# """
# Collects path scores for index positions on each chromosome.
# """
# def get_scores(
#     indices,
#     chromosomes,
#     score_type,
#     filter_orders,
#     full_file_path=os.path.join(datadir, 'scores/path_scores.order_14_resolution.h5'),
#     factor_file_path=os.path.join(datadir, 'scores/factor_scores.order_14_resolution.h5')
# ):
#     start = time.time()
#      # Preallocate scores array
#     scores = np.zeros(len(indices), dtype=np.float32)
#     print(f"scores for {len(indices)} indices")
#     # load in h5 file and read in the scores for each position
#     file_path = full_file_path if score_type == "full" else factor_file_path
#     with h5py.File(file_path, 'r') as f:
#         unique_chromosomes = np.unique(chromosomes)
#         for chromosome in unique_chromosomes:
#             csstart = time.time()
#             # chromosome_indices = indices[chromosomes == chromosome]
#             mask = chromosomes == chromosome
#             chromosome_indices = indices[mask]
#             print("chr indices", chromosome, len(indices[mask]))
#             sort_order = np.argsort(chromosome_indices)
#             unsort_order = np.argsort(sort_order)
#             print(f"Time to sort chromosome: {time.time() - csstart}", flush=True)

#             sstart = time.time()
#             # scores for each segment
#             chromosome_scores = f[chromosome][chromosome_indices[sort_order]][unsort_order]
#             print(f"Time to get chromosome scores: {time.time() - sstart}", flush=True)
            
#             # average scores for factors at each order filtered upon if not path scores
#             if score_type == "factor" and len(filter_orders) > 0:
#                 astart = time.time()
#                 orders_adjusted = np.array(filter_orders) - 4
#                 chromosome_scores = np.mean(chromosome_scores[:, orders_adjusted], axis=1)
#                 print(f"Time to average chromosome scores: {time.time() - astart}", flush=True)

#             # scores.extend(chromosome_scores)
#             scores[mask] = chromosome_scores

#     print(f"Time to get scores: {time.time() - start}", flush=True)
#     return np.array(scores)


"""
Refactors our top paths into a JSON object for return.
"""
def generate_return_data(chromosome_indices, top_indices, top_chromosomes, top_scores, dehydrated_paths, top_seen_counts):
    # return the top csns
    start = time.time()
    return_list = []
    # TODO: make this faster?
    file_path = os.path.join(datadir, 'scores/factor_scores.order_14_resolution.h5')
    with h5py.File(file_path, 'r') as f:
        for index, chrom, score, path, seen_count in zip(top_indices, top_chromosomes, top_scores, dehydrated_paths, top_seen_counts):
            scores = f[chrom][index]
            return_list.append({
                "i": int(index),
                "chromosome": str(chrom),
                "score": float(score),
                "path": path.tolist(),
                "factors": scores.tolist(),
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
    
    filters = request_info['filters']  # [{order, index, dataset_name}, ...]
    region = request_info['region']  # {order, chromosome, index}
    score_type = request_info['scoreType']  # "full", "factor"
    diversity = request_info['diversity']  # true | false
    N = request_info['N']  # number of paths to return

    
    # get the indices for the filters of each chromosome
    # if we have regions, only get the indices for chromomsomes present in the regions
    chromosome_indices = get_filtered_indices_and_scores(filters, region)
    # filter down the indices to the regions, converts indices to order of regions and checks
    if region is not None and len(filters) > 0:
        chromosome_indices = filter_region(chromosome_indices, region)

    filter_orders = [f['order'] for f in filters]
    max_order = max(filter_orders) if (diversity and len(filter_orders) > 0) else -1
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


"""
Endpoint for previewing the number of regions a filter would provide at each order if selected
"""
@csns_bp.route("preview_filter", methods=['POST'])
def preview_filter():
    request_info = request.get_json()
    
    filters = request_info['filters']  # [{order, index, dataset_name}, ...]
    new_filter = request_info['newFilter']  # [{order, index}, ...]
    region = request_info['region']  # {order, chromosome, index}

    # print((filters, new_filter, region), flush=True)

    
    # get the indices for the filters of each chromosome
    # if we have regions, only get the indices for chromomsomes present in the regions
    chromosome_indices = get_filtered_indices_and_scores(filters, region)
    print(chromosome_indices, flush=True)
    # # filter down the indices to the regions, converts indices to order of regions and checks
    # if region is not None and len(filters) > 0:
    #     chromosome_indices = filter_region(chromosome_indices, region)

    # filter_orders = [f['order'] for f in filters]
    # max_order = max(filter_orders) if (diversity and len(filter_orders) > 0) else -1
    # # find the top paths based on full or factor score
    # # print("sorting by score", flush=True)
    # top_indices, top_chromosomes, top_scores, top_seen_counts = sort_by_score(chromosome_indices, score_type, N, filter_orders, order_diversity=max_order)      
    
    # # get the dehydrated paths for each index 
    # # print("getting dehydrated paths", flush=True)
    # dehydrated_paths = get_dehydrated_paths(top_indices, top_chromosomes)

    # # generate the return list
    # # print("generating return list", flush=True)
    # return_list = generate_return_data(chromosome_indices, top_indices, top_chromosomes, top_scores, dehydrated_paths, top_seen_counts)

    return #return_list