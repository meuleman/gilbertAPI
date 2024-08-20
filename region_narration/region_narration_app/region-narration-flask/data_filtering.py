from flask import Blueprint, request
import numpy as np
import requests
import json
import h5py
import os
import time
from filtering import fetch_indices_scores

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

datadir="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/data_filtering/20240819"

def hilbert_pos_to_order(pos, frm, to):
    return np.floor(pos / pow(4, frm - to)).astype(int)


"""
For a given chromosome group, find the common indices between all filters.
"""


#TODO: Need to convert positions from one order to another (smaller order to larger order makes the most sense) before intersecting


def find_intersection(group):
    indices = group[0]["indices"]
    indices_sort_order = np.argsort(indices)
    indices = indices[indices_sort_order]
    combined_scores = group[0]["scores"][indices_sort_order]
    
    # add other groups to the combined indices and scores
    for g in group[1:]:
        g_indices = g["indices"]
        g_sort_order = np.argsort(g_indices)
        g_indices = g_indices[g_sort_order]
        g_factor_scores = g["scores"][g_sort_order]
        
        # find overlap in indices
        common_mask = np.isin(indices, g_indices)
        indices = indices[common_mask]
        combined_scores = combined_scores[common_mask]
       
       # combine factor scores
        g_mask = np.isin(g_indices, indices)
        combined_scores += g_factor_scores[g_mask]
    
    # average combined factor scores
    combined_scores = combined_scores / len(group)
    return {"indices": indices, "scores": combined_scores}


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
                os.path.join(datadir, 'indices.scores.order_14_resolution.h5')
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


"""
Endpoint for fetching the native order resolution region indices and data for any filter at any order
"""
@csns_bp.route("data_filtering", methods=['POST'])
def data_filtering():
    request_info = request.get_json()
    
    filters = request_info['filters']  # [{order, index, dataset_name}, ...]

    # find the intersection of filters
    chromosome_indices = get_filtered_indices_scores(filters)
    return