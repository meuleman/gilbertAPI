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

ORDER_SEGMENT_COUNTS = {
    4: 2957,
    5: 11792,
    6: 47136,
    7: 188505,
    8: 753983,
    9: 3015901,
    10: 12063568,
    11: 48254229,
    12: 193016875,
    13: 772067468,
    14: 3088269856
}

datadir="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/data_filtering/20240819"

"""
Function to preload all indices to be used for filter previews up through the max order
"""
def preload_preview_indices(max_order=12, index_path=os.path.join(datadir, 'indices.scores.native_resolution.h5')):
    preload_time = time.time()
    preloaded_indices = {}
    
    with h5py.File(index_path, 'r') as f:

        # find all the datasets less than or equal to the max order
        datasets = list(f.keys())
        datasets = [d for d in datasets if int(d.split('.')[0]) <= max_order]
        
        # preload indices for each dataset
        for d in datasets:

            # create chromosome independent key
            order, chrm, dataset_name, factor_index = d.split('.')
            key = f"{order}.{dataset_name}.{factor_index}"

            # all factor instances across chromosomes are saved together in a list
            if key not in preloaded_indices:
                preloaded_indices[key] = []
            
            # add to factorXorder list
            preloaded_indices[key].append({
                "chromosome": chrm,
                "order": int(order), 
                "index": int(factor_index), 
                "dataset_name": dataset_name, 
                "indices": f[d][:, 0].astype("int32"), 
                "scores": f[d][:, 1].astype("float32")
            })

    # save the order through which indices have been preloaded
    preloaded_indices['max_order'] = max_order
            
    print(f"Time to preload indices for filter previews: {time.time() - preload_time}", flush=True)
    return preloaded_indices

# preload all indices for filter previews
PRELOADED_INDICES = preload_preview_indices()



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
    combined_scores = group[0]["scores"][indices_sort_order] if "scores" in group[0] else None
    score_counter = sum([1 for g in group if "scores" in g])
    current_order = group[0]['order']

    # add other groups to the combined indices and scores
    for g in group[1:]:
        g_indices = g["indices"]
        g_sort_order = np.argsort(g_indices)
        g_indices = g_indices[g_sort_order]
        g_scores = g["scores"][g_sort_order] if "scores" in g else None

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
        combined_scores = combined_scores[common_mask] if combined_scores is not None else None

        # combine scores
        g_mask = np.isin(g_indices_base, indices_g_order)
        g_indices_base = g_indices_base[g_mask]
        g_scores = g_scores[g_mask] if g_scores is not None else None

        # count how many times each base position is found with new order indices
        g_indices_base_unique, g_indices_base_counts = np.unique(g_indices_base, return_counts=True)
        # expand the combined scores to match the position occurrences of the new order
        if np.count_nonzero(g_indices_base_unique != indices_g_order) > 0:
            raise ValueError(f"Inconsistent mapping between order {current_order} and order {g['order']} for {g['chromosome']} when calculating intersection")
        combined_scores = np.repeat(combined_scores, g_indices_base_counts) if combined_scores is not None else None

        indices = g_indices[g_mask]
        if ((combined_scores is not None) & (g_scores is not None)):
            combined_scores += g_scores
        elif (combined_scores is not None):
            combined_scores = combined_scores
        else:
            combined_scores = g_scores

        # update current order
        current_order = g['order']

    # average combined factor scores
    combined_scores = combined_scores / score_counter
    return {"indices": indices, "scores": combined_scores, 'order': current_order}


def get_filtered_indices_scores(filters, force_h5=True):
    filter_indices = []
    start = time.time()
    
    # parse each filter and collect indices
    for f in filters:
        order = f['order']
        index = f['index']
        dataset_name = f['dataset_name']

        # create chromosome independent key
        key = f"{order}.{dataset_name}.{index}"

        # if not force_h5, fetch preloaded indices
        if ((not force_h5) & (int(order) <= PRELOADED_INDICES['max_order'])):
            filter_indices.extend(PRELOADED_INDICES[key]) if key in PRELOADED_INDICES else None

        else:
            # if force_h5, fetch indices from h5 file for each chromosome
            for chromosome in chromosomes:

                # load indices and scores
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

    # print(f"Time to fetch indices: {time.time() - start}", flush=True)

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
    
    # print(f"Time to find intersection: {time.time() - start}", flush=True)
    return reduced_chromosome_groups


"""
Using the combined factor scores for index positions across chromosomes, sorts by score
and returns the indices, chromosomes, and scores for the top N segments across the genome.
"""
def sort_by_score(chromosome_indices, N):
    start = time.time()

    # initialize arrays
    indices = np.array([], dtype="int32")
    scores = np.array([], dtype="float32")
    chromosomes = np.array([])

    # collect indices and scores from each chromosome
    for chromosome in chromosome_indices.keys():
        chrm_indices = chromosome_indices[chromosome]["indices"]
        chrm_scores = chromosome_indices[chromosome]['scores']
        indices = np.append(indices, chrm_indices)
        scores = np.append(scores, chrm_scores)
        chromosomes = np.append(chromosomes, np.array([chromosome] * len(chrm_indices)))

    # count number of segments
    segment_count = len(indices)

    # sort by score
    sort_order = np.argsort(-scores)[:N]
    top_scores = scores[sort_order]
    top_indices = indices[sort_order]
    top_chromosomes = chromosomes[sort_order]

    print(f"Time to sort by score: {time.time() - start}", flush=True)

    return top_indices, top_chromosomes, top_scores, segment_count


"""
Refactors our top segments into a JSON object for return.
"""
def generate_return_data(top_indices, top_chromosomes, top_scores, segment_count, chromosome_indices):
    start = time.time()

    # # return the filtered indices
    # return_indices = []
    # for chr in chromosome_indices.keys():
    #     return_indices.append({
    #         "chromosome": chr,
    #         "indices": chromosome_indices[chr]['indices'].tolist(),
    #         "scores": chromosome_indices[chr]['scores'].tolist(),
    #         "order": chromosome_indices[chr]['order']
    #     })


    # return the top segments
    return_segments = []
    for index, chrm, score in zip(top_indices, top_chromosomes, top_scores):
        return_segments.append({'chromosome': chrm, 'index': int(index), 'score': float(score)})

    print(f"Time to generate return data: {time.time() - start}", flush=True)
    return json.dumps({
        "filtered_segments": return_segments,
        "order": list(chromosome_indices.values())[0]['order'] if len(chromosome_indices) > 0 else None,
        "segment_count": segment_count
    })



"""
Finds the regional index overlap between established and proposed filters and returns fractional overlap value at each order
"""
def get_factor_filter_previews(regions, new_filters, orders=range(4, 15), normalize=True):
    # get the indices for the regions we have established for each chromosome
    established_start = time.time()
    set_region_order = None
    if len(regions) > 0:  # if no regions in place, skip
        set_chromosome_indices = regions_to_chrm_groups(regions)
        # find the chromosomes for set regions
        set_chromosomes = np.array(list(set_chromosome_indices.keys()))
        # to count the total number of set regions at each new order later
        set_indices = [set_chromosome_indices[chromosome]['indices'] for chromosome in set_chromosomes]
        # order resolution for the set regions
        set_region_order = set_chromosome_indices[set_chromosomes[0]]['order'] if len(set_chromosomes) > 0 else None
    print(f"Time to get established filter regions: {time.time() - established_start}", flush=True)

    # loop through each new filter and find overlap at each order
    previews = []
    for new_filter in new_filters:
        new_filter_start = time.time()
        # loop through each order and find the overlap between set filters and new filter
        preview = {}
        for order in orders:
            new_filter_at_order = new_filter.copy()
            new_filter_at_order['order'] = order
            # get the regional indices for the new filter
            preview_chromosome_indices = get_filtered_indices_scores([new_filter_at_order], force_h5=False)
            # find the chromosomes new filter has regions for
            preview_chromosomes = np.array(list(preview_chromosome_indices.keys()))

            # for order, find the common chromosomes for set regions and new filters
            commom_chromosomes = set_chromosomes[np.isin(set_chromosomes, preview_chromosomes)] if len(regions) > 0 else preview_chromosomes

            # group indices between filters by chromosome and find intersection counts
            common_indices = []
            for chromosome in commom_chromosomes:
                if len(regions) > 0:  # find intersection if there are set regions
                    chromosome_group = [set_chromosome_indices[chromosome], preview_chromosome_indices[chromosome]]
                    chromosome_intersection = find_intersection(chromosome_group)
                    common_indices.append(chromosome_intersection['indices'])
                else:
                    common_indices.append(preview_chromosome_indices[chromosome]['indices'])
                
            # convert common indices to the preview order
            max_common_order = max(set_region_order, order) if set_region_order != None else order
            unique_order_common_indices = [
                np.unique(hilbert_pos_to_order(indices, max_common_order, order)) for indices in common_indices
            ]
            common_count = sum([len(indices) for indices in unique_order_common_indices])
            
            if len(regions) > 0:
                # find number of unique indices from set regions at order
                unique_order_set_indices = [
                    np.unique(hilbert_pos_to_order(indices, set_region_order, order)) for indices in set_indices
                ]
                set_count = sum([len(indices) for indices in unique_order_set_indices])
                # if the order for set regions is < the preview order, we need to adjust our count
                if set_region_order != None:
                    if set_region_order < max_common_order:
                        set_count = set_count * (4 ** (max_common_order - set_region_order))
                
            
            # if no set regions use full genome size as background, otherwise use regions
            bg_count = set_count if len(regions) > 0 else ORDER_SEGMENT_COUNTS[order]
            # print("ORDER", order, "COUNT", common_count, "BG", bg_count)
            if normalize:
                preview[order] = common_count / bg_count if bg_count > 0 else 0
            else:
                preview[order] = common_count
        
        filter_preview = {'preview': preview, 'dataset_name': new_filter['dataset_name'], 'index': new_filter['index']}
        previews.append(filter_preview)
        # print(f"Time to find preview for filter {new_filter['dataset_name']}, {new_filter['index']}: {time.time() - new_filter_start}", flush=True)
    
    return previews



# """
# Finds the regional index overlap between established and proposed filters and returns fractional overlap value at each order
# """
# def get_factor_filter_previews(filters, new_filters, orders=range(4, 15), normalize=True):
#     # get the indices for the filters we have established for each chromosome
#     established_start = time.time()
#     set_filter_order = None
#     if len(filters) > 0:  # if no filters in place, skip
#         set_chromosome_indices = get_filtered_indices_scores(filters, force_h5=False)
#         # find the chromosomes set filters have regions for
#         set_chromosomes = np.array(list(set_chromosome_indices.keys()))
#         # to count the total number of regions for set filters at each new order later
#         set_indices = [set_chromosome_indices[chromosome]['indices'] for chromosome in set_chromosomes]
#         # order resolution for the set filters
#         set_filter_order = set_chromosome_indices[set_chromosomes[0]]['order'] if len(set_chromosomes) > 0 else None
#     print(f"Time to get established filter regions: {time.time() - established_start}", flush=True)

#     # loop through each new filter and find overlap at each order
#     previews = []
#     for new_filter in new_filters:
#         new_filter_start = time.time()
#         # loop through each order and find the overlap between set filters and new filter
#         preview = {}
#         for order in orders:
#             new_filter_at_order = new_filter.copy()
#             new_filter_at_order['order'] = order
#             # get the regional indices for the new filter
#             preview_chromosome_indices = get_filtered_indices_scores([new_filter_at_order], force_h5=False)
#             # find the chromosomes new filter has regions for
#             preview_chromosomes = np.array(list(preview_chromosome_indices.keys()))

#             # for order, find the common chromosomes for set and new filters
#             commom_chromosomes = set_chromosomes[np.isin(set_chromosomes, preview_chromosomes)] if len(filters) > 0 else preview_chromosomes

#             # group indices between filters by chromosome and find intersection counts
#             common_indices = []
#             for chromosome in commom_chromosomes:
#                 if len(filters) > 0:  # find intersection if there are set filters
#                     chromosome_group = [set_chromosome_indices[chromosome], preview_chromosome_indices[chromosome]]
#                     chromosome_intersection = find_intersection(chromosome_group)
#                     common_indices.append(chromosome_intersection['indices'])
#                 else:
#                     common_indices.append(preview_chromosome_indices[chromosome]['indices'])
                
#             # convert common indices to the preview order
#             max_common_order = max(set_filter_order, order) if set_filter_order != None else order
#             unique_order_common_indices = [
#                 np.unique(hilbert_pos_to_order(indices, max_common_order, order)) for indices in common_indices
#             ]
#             common_count = sum([len(indices) for indices in unique_order_common_indices])
            
#             if len(filters) > 0:
#                 # find number of unique indices from set filters at order
#                 unique_order_set_indices = [
#                     np.unique(hilbert_pos_to_order(indices, set_filter_order, order)) for indices in set_indices
#                 ]
#                 set_count = sum([len(indices) for indices in unique_order_set_indices])
#                 # if the order for set filters is < the preview order, we need to adjust our count
#                 if set_filter_order != None:
#                     if set_filter_order < max_common_order:
#                         set_count = set_count * (4 ** (max_common_order - set_filter_order))
                
            
#             # if no set filters, use full genome size as background, otherwise use regions found through filters
#             bg_count = set_count if len(filters) > 0 else ORDER_SEGMENT_COUNTS[order]
#             # print("ORDER", order, "COUNT", common_count, "BG", bg_count)
#             if normalize:
#                 preview[order] = common_count / bg_count if bg_count > 0 else 0
#             else:
#                 preview[order] = common_count
        
#         filter_preview = {'preview': preview, 'dataset_name': new_filter['dataset_name'], 'index': new_filter['index']}
#         previews.append(filter_preview)
#         # print(f"Time to find preview for filter {new_filter['dataset_name']}, {new_filter['index']}: {time.time() - new_filter_start}", flush=True)
    
#     return previews


"""
Converts a list of regions to chromosome groups
"""
def regions_to_chrm_groups(regions):
    region_chromosome_indices = {}
    for region in regions:
        if region['chromosome'] not in region_chromosome_indices:
            region_chromosome_indices[region['chromosome']] = {
                "chromosome": region['chromosome'], 
                "order": region['order'],
                "indices": np.array([]).astype(int),
            }
        region_chromosome_indices[region['chromosome']]['indices'] = np.append(
            region_chromosome_indices[region['chromosome']]['indices'], region['index']
        )

    return region_chromosome_indices


"""
Takes in a dictionary of chromosome indices and filters them by regions.
This assumes all input regions exist at the same order and all filter regions
exist at the same order.
"""
def filter_region(chromosome_indices, regions):

    # convert regions to chromosome groups
    region_chromosome_indices = regions_to_chrm_groups(regions)

    # group by chromosome
    chromosome_groups = {}
    for chrm in chromosome_indices:
        # only keep if regions on this chromosome
        if chrm in region_chromosome_indices:
            chromosome_group = [chromosome_indices[chrm], region_chromosome_indices[chrm]]
            chromosome_groups[chrm] = chromosome_group


    # reduce chromosome groups to common indices
    reduced_chromosome_groups = {}
    for chromosome, group in chromosome_groups.items():
        print(("GROUP", group), flush=True)
        group_intersection = find_intersection(group)
        if len(group_intersection["indices"]) > 0:
            reduced_chromosome_groups[chromosome] = group_intersection

    return reduced_chromosome_groups



"""
Endpoint for fetching the native order resolution region indices and data for any filter at any order
"""
@data_filtering_bp.route("data_filtering", methods=['POST'])
def data_filtering():
    request_info = request.get_json()
    
    filters = request_info['filters']  # [{order, index, dataset_name}, ...]
    regions = request_info['regions'] if 'regions' in request_info else None  # [{chromosome, order, index}, ...]
    N = request_info['N'] if 'N' in request_info else None  # number of segments to return

    # find the intersection of filters
    chromosome_indices = get_filtered_indices_scores(filters)

    # filter the indices by regions
    chromosome_indices = filter_region(chromosome_indices, regions) if regions is not None else chromosome_indices

    # find the top segments based on combined factor score
    top_indices, top_chromosomes, top_scores, segment_count = sort_by_score(chromosome_indices, N)

    # generate the return list
    return_list = generate_return_data(top_indices, top_chromosomes, top_scores, segment_count, chromosome_indices)
    
    return return_list


"""
Endpoint for previewing the number of segments a filter would provide at each order if selected
"""
@data_filtering_bp.route("preview_filter", methods=['POST'])
def preview_filter():
    request_info = request.get_json()
    print("REQUEST:", request.method)
    
    filters = request_info['filters']  # [{order, index, dataset_name}, ...]
    regions = request_info['regions'] if 'regions' in request_info else [] # [{chromosome, order, index}, ...]
    new_filters = request_info['newFilters']  # [{dataset_name, index}, ...]
    orders = request_info['orders'] if 'orders' in request_info else range(4, 15)  # [4, 5, 6, ...]
    normalize = request_info['normalize'] if 'normalize' in request_info else True
    normalize = normalize.lower() == "true" if type(normalize) == str else normalize

    print(regions, flush=True)

    preview_start = time.time()
    previews = get_factor_filter_previews(regions, new_filters, orders, normalize)
    print("total time to find overlap of all factors", time.time() - preview_start, flush=True)

    return_object = json.dumps({"previews": previews})
    return return_object