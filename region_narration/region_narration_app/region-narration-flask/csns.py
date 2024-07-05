from flask import Blueprint, Flask, request, abort, Response
import numpy as np
import requests
import os
import json
from functools import reduce
import h5py

# Create a Blueprint
csns_bp = Blueprint('csns_bp', __name__)

baseS3URLPrefix = "https://altius-gilbert.s3.us-west-2.amazonaws.com"
chromosomes = np.append([f"chr{c}" for c in range(1, 23)], ['chrX', 'chrY'])



# # load path scores from each chromosome
# scores_dir = '/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/precomputed_csn/path_scores'
scores_lookup = {}
# for chromosome in chromosomes:
#     chromosome_score_path = os.path.join(scores_dir, f"path_scores.{chromosome}.order_14_resolution.float32.bytes")
#     with open(chromosome_score_path, "rb") as f:
#         bytes_data = f.read()
#         f.close()
#     scores_lookup[chromosome] = np.frombuffer(bytes_data, dtype="float32")


def fetch_indices(order, chromosome, factor_index, dataset_name):
    url = f"{baseS3URLPrefix}/${order}.${chromosome}.${dataset_name}.${factor_index}.order_14_resolution.indices.int32.bytes"
    bytes_arr = requests.get(url)
    indices = np.frombuffer(bytes_arr, dtype="int32")
    return indices

def find_intersection(group):
    common_indices = group[0]["indices"]
    for g in group[1:]:
        common_indices = common_indices[np.isin(common_indices, g["indices"])]
    return common_indices

def get_filtered_indices(filters):
    filter_indices = []
    # parse each filter and collect indices
    for f in filters:
        order = f['order']
        index = f['index']
        dataset_name = f['dataset_name']
        # fetch indices from each chromosome
        for chromosome in chromosomes:
            indices = fetch_indices(order, chromosome, index, dataset_name)
            filter_indices.append({
                "chromosome": chromosome, 
                "order": order, 
                "index": index, 
                "dataset_name": dataset_name, 
                "indices": indices
            })

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
    
    return reduced_chromosome_groups



def sort_by_score(chromosome_indices, score_type, N):
    indices = np.array([])
    chromosomes = []
    scores = []

    for chromosome in chromosome_indices.keys():
        chrm_indices = chromosome_indices[chromosome]
        chrm_array = np.array([chromosome] * len(chrm_indices))
        
        if score_type == "full":
            chrm_scores = scores_lookup[chromosome][chrm_indices]

        indices = np.append(indices, chrm_indices)
        chromosomes = np.append(chromosomes, chrm_array)
        scores = np.append(scores, chrm_scores)

    # sort
    sort_order = np.argsort(-scores)[:N]
    top_scores = scores[sort_order]
    top_indices = indices[sort_order]
    top_chromosomes = chromosomes[sort_order]

    return top_indices, top_chromosomes, top_scores


def get_dehydrated_paths(
    top_indices,
    top_chromosomes,
    file_path = '/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/precomputed_csn/paths/path_factor_indices.order_14_resolution.h5'
):
    dehydrated_paths = []
    # load in h5 file and read in the exact dehydrated paths we need
    with h5py.File(file_path, 'r') as f:
        for chromosome, index in zip(top_chromosomes, top_indices):
            dehydrated_paths.append(f[chromosome][index])
    
    return dehydrated_paths


def generate_return_list(top_indices, top_chromosomes, top_scores, dehydrated_paths):
    return_list = []
    for index, chrom, score, path in zip(top_indices, top_chromosomes, top_scores, dehydrated_paths):
        return_list.append({
            "index": index,
            "chromosome": chrom,
            "score": score,
            "path": path.tolist()
        })

    return json.dumps(return_list)


@csns_bp.route("top_paths")
def top_paths():
    request_info = request.get_json()
    
    filters = request_info['filters']  # [{order, index}, ...]
    regions = request_info['regions']  # [{order, chromosome, index},...]
    score_type = request_info['scoreType']  # "full", "factor"
    diversity = request_info['diversity']  # true | false
    N = request_info['N']  # number of paths to return

    # get the indices for the filters of each chromosome
    chromosome_indices = get_filtered_indices(filters)

    #TODO: filter by regions that are passed in

    # find the top paths based on full or factor score
    top_indices, top_chromosomes, top_scores = sort_by_score(chromosome_indices, score_type, N)        
    
    # get the dehydrated paths for each index 
    dehydrated_paths = get_dehydrated_paths(top_indices, top_chromosomes)

    # generate the return list
    return_list = generate_return_list(top_indices, top_chromosomes, top_scores, dehydrated_paths)

    return return_list