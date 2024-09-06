from flask import Blueprint, request
import numpy as np
import requests
import json
import h5py
import os
import time
from filtering import fetch_indices_scores
# import tabix
from subprocess import Popen, PIPE

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
genedir="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/genes/20240906/"
gene_path = os.path.join(
    genedir, 
    "mart_export.protein_coding.basal_plus_extension.5kb_upstream.1kb_downstream.1000kb_distal.CT20240906.bed.bgz"
)
# columns for gene_path
# [
# 'Chromosome', 'Left distal', 'Right distal', 'Gene name', 
# 'Gene stable ID', 'Strand', 'Gene start (bp)', 'Gene end (bp)', 
# 'Basal upstream', 'Basal downstream', 'Gene description'
# ]

def hilbert_pos_to_order(pos, frm, to):
    return np.floor(pos / pow(4, frm - to)).astype(int)


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
        scores = get_single_region_scores(chromosome, start, end)
        print("len indices", len(indices), "len scores", len(scores))
        # no need to deal with other chromosomes
        return {
            chromosome: {"indices": indices, "path_scores": scores}
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
            indices_scores = fetch_indices_scores(
                order, 
                chromosome, 
                index, 
                dataset_name, 
                os.path.join(datadir, 'indices/indices.scores.order_14_resolution.h5')
            )
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
        group_intersection = find_intersection(group)
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
def sort_by_score(chromosome_indices, score_type, N, order_diversity=-1):
    indices = np.array([], dtype="int32")
    factor_scores = np.array([], dtype="float32")
    path_scores = np.array([], dtype="float32")
    # score_key = "factor_scores" if score_type == 'factor' else "path_scores"
    chromosomes = np.array([])
    seen_counts = np.array([], dtype="int32")

    start = time.time()
    for chromosome in chromosome_indices.keys():
        chrm_indices = chromosome_indices[chromosome]["indices"]
        if score_type == "factor":
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
            if score_type == "factor":
                chrm_factor_scores = chrm_factor_scores[chrm_ops_indices]
            chrm_path_scores = chrm_path_scores[chrm_ops_indices]
        else:
            chrm_ops_counts = np.ones(len(chrm_indices), dtype="int32")
        
        indices = np.append(indices, chrm_indices)
        if score_type == "factor":
            factor_scores = np.append(factor_scores, chrm_factor_scores)
        path_scores = np.append(path_scores, chrm_path_scores)
        chromosomes = np.append(chromosomes, np.array([chromosome] * len(chrm_indices)))
        seen_counts = np.append(seen_counts, chrm_ops_counts)

    print(f"Time to find diversity: {time.time() - start}", flush=True)

    start = time.time()
    # sort
    if score_type == "factor":
        sort_order = np.argsort(-factor_scores)[:N]
    else:
        sort_order = np.argsort(-path_scores)[:N]
    # We always return path scores, its just the sorting determined by the score type
    top_scores = path_scores[sort_order]
    top_indices = indices[sort_order]
    top_chromosomes = chromosomes[sort_order]
    top_seen_counts = seen_counts[sort_order]

    print(f"Time to sort by score: {time.time() - start}", flush=True)

    return top_indices, top_chromosomes, top_scores, top_seen_counts


"""
Collects and returns the factor indices for each of the top paths.
"""
def get_dehydrated_paths(
    top_indices=[],
    top_chromosomes=[],
    regions=[],
    file_path=os.path.join(datadir, 'paths/full_path_factors.order_14_resolution.h5')
):
    start = time.time()

    if len(regions) > 0:
        # load in h5 file and read in the exact dehydrated paths we need for each region
        with h5py.File(file_path, 'r') as f:
            for region in regions:
                paths = []
                for pos in region['top_positions']:
                    paths.append(f[region['chromosome']][pos].tolist())
                region['dehydrated_paths'] = paths

        print(f"Time to get dehydrated paths: {time.time() - start}", flush=True)
        return regions
    
    elif (len(top_indices) > 0) and (len(top_chromosomes) > 0):
        dehydrated_paths = []
        # load in h5 file and read in the exact dehydrated paths we need
        with h5py.File(file_path, 'r') as f:
            for chromosome, index in zip(top_chromosomes, top_indices):
                dehydrated_paths.append(f[chromosome][index])
        
        print(f"Time to get dehydrated paths: {time.time() - start}", flush=True)
        return dehydrated_paths
    
    else:
        raise ValueError("No indices or chromosomes provided for dehydrated paths")


"""
Collects path scores for index positions on each chromosome.
"""
def get_scores(
    indices,
    chromosomes,
    score_type,
    filter_orders,
    full_file_path=os.path.join(datadir, 'scores/path_scores.order_14_resolution.h5'),
    factor_file_path=os.path.join(datadir, 'scores/factor_scores.order_14_resolution.h5')
):
    start = time.time()
    # Preallocate scores array
    scores = np.zeros(len(indices), dtype=np.float32)
    print(f"scores for {len(indices)} indices")
    # load in h5 file and read in the scores for each position
    file_path = full_file_path if score_type == "full" else factor_file_path
    with h5py.File(file_path, 'r') as f:
        unique_chromosomes = np.unique(chromosomes)
        for chromosome in unique_chromosomes:
            csstart = time.time()
            # chromosome_indices = indices[chromosomes == chromosome]
            mask = chromosomes == chromosome
            chromosome_indices = indices[mask]
            print("chr indices", chromosome, len(indices[mask]))
            sort_order = np.argsort(chromosome_indices)
            unsort_order = np.argsort(sort_order)
            print(f"Time to sort chromosome: {time.time() - csstart}", flush=True)

            sstart = time.time()
            # scores for each segment
            chromosome_scores = f[chromosome][chromosome_indices[sort_order]][unsort_order]
            print(f"Time to get chromosome scores: {time.time() - sstart}", flush=True)
            
            # average scores for factors at each order filtered upon if not path scores
            if score_type == "factor" and len(filter_orders) > 0:
                astart = time.time()
                orders_adjusted = np.array(filter_orders) - 4
                chromosome_scores = np.mean(chromosome_scores[:, orders_adjusted], axis=1)
                print(f"Time to average chromosome scores: {time.time() - astart}", flush=True)

            # scores.extend(chromosome_scores)
            scores[mask] = chromosome_scores

    print(f"Time to get scores: {time.time() - start}", flush=True)
    return np.array(scores)


"""
Collects path scores for a single region on a chromosome.
"""
def get_single_region_scores(
    chromosome,
    start,
    end,
    file_path=os.path.join(datadir, 'scores/path_scores.order_14_resolution.h5'),
):
    score_t = time.time()
    print(f"scores for {(end - start)} indices")
    # load in h5 file and read in the scores for each position
    with h5py.File(file_path, 'r') as f:
        # scores for each segment
        scores = f[chromosome][start:end]

    print(f"Time to get scores: {time.time() - score_t}", flush=True)
    return scores


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
    #         "indices": chromosome_indices[chr]["indices"][0:1000000].tolist()
    #     })

    print(f"Time to generate return data: {time.time() - start}", flush=True)
    return json.dumps({
        # "filtered_indices": return_indices,
        "csns": return_list
    })

"""
Finds the regional index overlap between established and proposed filters and returns fractional overlap value at each order
"""
def get_factor_filter_previews(filters, new_filter, region, orders=range(4, 15)):
    # get the indices for the filters we have established for each chromosome
    # if we have regions, only get the indices for chromomsomes present in the regions
    if len(filters) > 0:  # if no filters in place, skip
        set_chromosome_indices = get_filtered_indices_and_scores(filters, region)
        # find the chromosomes set filters have regions for
        set_chromosomes = np.array(list(set_chromosome_indices.keys()))
        # count the total number of regions for set filters
        set_index_count = sum(len(set_chromosome_indices[chromosome]['indices']) for chromosome in set_chromosomes)

    # loop through each order and find the overlap between set filters and new filter
    previews = {}
    for order in orders:
        new_filter_at_order = new_filter.copy()
        new_filter_at_order['order'] = order
        # get the regional indices for the new filter
        preview_chromosome_indices = get_filtered_indices_and_scores([new_filter_at_order], region)
        # find the chromosomes new filter has regions for
        preview_chromosomes = np.array(list(preview_chromosome_indices.keys()))

        # for order, find the common chromosomes for set and new filters
        commom_chromosomes = set_chromosomes[np.isin(set_chromosomes, preview_chromosomes)] if len(filters) > 0 else preview_chromosomes

        # group indices between filters by chromosome and find intersection counts
        common_index_count = 0
        for chromosome in commom_chromosomes:
            if len(filters) > 0:  # find intersection if there are set filters
                chromosome_group = [set_chromosome_indices[chromosome], preview_chromosome_indices[chromosome]]
                chromosome_intersection = find_intersection(chromosome_group)
                common_index_count += len(chromosome_intersection['indices'])
            else:
                common_index_count += len(preview_chromosome_indices[chromosome]['indices'])

        # if no set filters, use full genome size as background, otherwise use regions found through filters
        bg_count = set_index_count if len(filters) > 0 else sum(CHR_SIZES_UCSC.values())
        print("ORDER", order, "COUNT", common_index_count, "BG", bg_count)
        previews[order] = common_index_count / bg_count

    return previews


"""
Collect the top N scores and positions for each region
"""
def get_top_scores_for_regions(regions, N):
    start_time = time.time()
    print(f"Finding top scores for {len(regions)} regions")

    for region in regions:
        top_indices = np.argsort(-region['scores'])[:N]
        top_scores = region['scores'][top_indices]
        top_positions = np.arange(region['start'], region['end'])[top_indices]
        region['top_scores'] = top_scores
        region['top_positions'] = top_positions

    print(f"Time to find top scores: {time.time() - start_time}", flush=True)

    return regions


"""
Collects path scores for each region set. 
"""
def get_scores_for_regions(
    regions,
    file_path=os.path.join(datadir, 'scores/path_scores.order_14_resolution.h5'),
):
    start_time = time.time()
    print(f"Collecting scores for {len(regions)} regions")

    # load in h5 file and read in the scores for each position
    with h5py.File(file_path, 'r') as f:
        for region in regions:
            region['scores'] = f[region['chromosome']][region['start']:region['end']]

    print(f"Time to collect scores: {time.time() - start_time}", flush=True)

    return regions


"""
Create return object for regions
"""
def generate_region_return(regions):
    for region in regions:
        del region['scores']
        region['top_scores'] = region['top_scores'].tolist()
        region['top_positions'] = region['top_positions'].tolist()
        region['dehydrated_paths'] = region['dehydrated_paths']

    return json.dumps({"regions": regions})


"""
Queries a tabix file for provided genomic coordinates (chr:start-end)
Return all the overlapping hits
"""
def tabix_query(query, filename):
    # Use subprocess with tabix and return hits
    process = Popen(['tabix', '-f', filename, query], stdout=PIPE)
    
    hits = []
    for line in process.stdout:
        hits.append(str(np.char.decode(line)).strip().split('\t'))
    return hits


"""
Collects the genes associated with each path in each region.
"""
def get_genes(regions, filepath=gene_path):
    start = time.time()
    # for each region
    for region in regions:
        region_genes = []
        for pos in region['top_positions']:
            pos_coords = f"{region['chromosome']}:{pos}-{pos}"
            genes = tabix_query(pos_coords, filepath)
            region_genes.append({
                'i': pos,
                'genes': [
                        {
                            'chromosome': g[0], 
                            'domain_start': int(g[1]),
                            'domain_end': int(g[2]),
                            'name': g[3],
                            'stable_ID': g[4],
                            'strand': int(g[5]),
                            'start': int(g[6]),
                            'end': int(g[7]),
                            'description': g[10],
                            'in_gene': ((pos >= int(g[6])) & (pos <= int(g[7])))
                        } for g in genes
                ]
            })
        region['genes'] = region_genes

    print(f"Time to collect gene information for regions: {time.time() - start}", flush=True)
    return regions


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
    top_indices, top_chromosomes, top_scores, top_seen_counts = sort_by_score(chromosome_indices, score_type, N, order_diversity=max_order)      
    
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

    print(filters, new_filter)

    previews = get_factor_filter_previews(filters, new_filter, region)

    return_object = json.dumps({"preview_fractions": previews})
    return return_object


"""
Endpoint for fetching the top paths for a number of regions. The top N paths for each region are returned.
"""
@csns_bp.route("top_paths_for_regions", methods=['POST'])
def native_order_resolution_regions():
    request_info = request.get_json()
    
    region_coordinates = request_info['regions']  # ['chr:start-end', ...]
    N = request_info['N'] if 'N' in request_info else 1 # number of paths to return per region (default 1)

    # convert region coordinates to regions object
    regions = [
        {
            'chromosome': r.split(':')[0],
            'start': int(r.split(':')[1].split('-')[0]),
            'end': int(r.split(':')[1].split('-')[1])
        } for r in region_coordinates
    ]

    # get the scores for each region
    regions = get_scores_for_regions(regions)

    # find the top N positions based on path score
    regions = get_top_scores_for_regions(regions, N)

    # get the genes for each position in each region
    regions = get_genes(regions)

    # get the dehydrated paths for each position in each region
    regions = get_dehydrated_paths(regions=regions)

    # generate the return object
    return_object = generate_region_return(regions)
    
    return return_object