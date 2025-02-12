from flask import Blueprint, request
import numpy as np
# import requests
import json
import h5py
import os
import time
from common import fetch_indices_scores, hilbert_pos_to_order, tabix_query
# import tabix
from constants import CHR_SIZES_UCSC, chromosomes, csn_dir, full_gene_path, partial_csn_dir

# Create a Blueprint
csns_bp = Blueprint('csns_bp', __name__)

# baseS3URLPrefix = "https://altius-gilbert.s3.us-west-2.amazonaws.com/20240703/csn_index_files"  # dataset removed from s3


def find_intersection(group):
    """
    For a given chromosome group, find the common indices between all filters.
    """
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


def get_filtered_indices_and_scores(filters, region):
    """
    Collects positional indices for each filter, groups them by chromosome, and finds the
    intersection between filters within each chromosome.
    """
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
                os.path.join(csn_dir, 'indices/indices.scores.order_14_resolution.h5')
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

def filter_region(chromosome_indices, region):
    """
    Takes in a dictionary of chromosome indices and filters them by regions.
    If more than one regions are present, we OR the regions together, giving all paths that flow through each region
    """
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

def sort_by_score(chromosome_indices, score_type, N, order_diversity=-1):
    """
    Collects the path scores for index positions across chromosome, sorts by score,
    and returns the indices, chromosomes, and scores for the top N paths across the genome.
    """
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


def get_dehydrated_paths(
    top_indices=[],
    top_chromosomes=[],
    regions=[],
    file_path=os.path.join(csn_dir, 'paths/full_path_factors.order_14_resolution.h5')
):
    """
    Collects and returns the factor indices for each of the top paths.
    """
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


def get_scores(
    indices,
    chromosomes,
    score_type,
    filter_orders,
    full_file_path=os.path.join(csn_dir, 'scores/path_scores.order_14_resolution.h5'),
    factor_file_path=os.path.join(csn_dir, 'scores/factor_scores.order_14_resolution.h5')
):
    """
    Collects path scores for index positions on each chromosome.
    """
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


def get_single_region_scores(
    chromosome,
    start,
    end,
    file_path=os.path.join(csn_dir, 'scores/path_scores.order_14_resolution.h5'),
):
    """
    Collects path scores for a single region on a chromosome.
    """
    score_t = time.time()
    print(f"scores for {(end - start)} indices")
    # load in h5 file and read in the scores for each position
    with h5py.File(file_path, 'r') as f:
        # scores for each segment
        scores = f[chromosome][start:end]

    print(f"Time to get scores: {time.time() - score_t}", flush=True)
    return scores


def generate_return_data(chromosome_indices, top_indices, top_chromosomes, top_scores, dehydrated_paths, top_seen_counts):
    """
    Refactors our top paths into a JSON object for return.
    """
    # return the top csns
    start = time.time()
    return_list = []
    # TODO: make this faster?
    file_path = os.path.join(csn_dir, 'scores/factor_scores.order_14_resolution.h5')
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

def get_factor_filter_previews(filters, new_filter, region, orders=range(4, 15)):
    """
    Finds the regional index overlap between established and proposed filters and returns fractional overlap value at each order
    """
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


def get_top_scores_for_regions(regions, N):
    """
    Collect the top N scores and positions for each region
    """
    start_time = time.time()
    print(f"Finding top scores for {len(regions)} regions")

    for region in regions:
        top_indices = np.argsort(-region['path_scores'])[:N]
        top_path_scores = region['path_scores'][top_indices]
        top_positions = np.arange(region['start'], region['end'])[top_indices]
        region['top_path_scores'] = top_path_scores
        region['top_positions'] = top_positions

    print(f"Time to find top scores: {time.time() - start_time}", flush=True)

    return regions


def get_scores_for_regions(
    regions,
    score_type,
    factor_score_path=os.path.join(csn_dir, 'scores/factor_scores.order_14_resolution.h5'),
    path_score_path=os.path.join(csn_dir, 'scores/path_scores.order_14_resolution.h5'),
    by_position=False
):
    """
    Collects path scores for each region set. 
    """
    start_time = time.time()
    print(f"Collecting scores for {len(regions)} regions")

    # load in h5 file and read in the scores for each position
    file_path, region_key = (path_score_path, "path_scores") if score_type == "path" else (factor_score_path, "factor_scores")
    with h5py.File(file_path, 'r') as f:
        for region in regions:
            # if collecting scores for top positions in region
            if by_position:
                region_key = "top_path_scores" if score_type == "path" else "top_factor_scores"
                pos_scores = []
                for pos in region['top_positions']:
                    pos_scores.append(f[region['chromosome']][pos])
                region[region_key] = np.array(pos_scores)
            # if collecting scores for the entire region
            else:
                region[region_key] = f[region['chromosome']][region['start']:region['end']]

    print(f"Time to collect scores: {time.time() - start_time}", flush=True)

    return regions


def generate_region_return(regions):
    """
    Create return object for regions
    """
    for region in regions:
        if 'path_scores' in region:
            del region['path_scores']
        if 'factor_scores' in region:
            del region['factor_scores']
        region['top_path_scores'] = region['top_path_scores'].tolist()
        region['top_factor_scores'] = region['top_factor_scores'].tolist()
        region['top_positions'] = region['top_positions'].tolist()

    return json.dumps({"regions": regions})


def get_genes(chrm, start, end):
    """
    Given chromosome, start, and end, finds the genes overlapping the region.
    """
    coords = f"{chrm}:{start}-{end}"
    genes = tabix_query(coords, full_gene_path)
    genes_arr = [
        {
            'chromosome': str(g[0]), 
            'domain_start': int(g[1]),
            'domain_end': int(g[2]),
            'name': str(g[3]),
            'stable_ID': str(g[4]),
            'strand': int(g[5]),
            'start': int(g[6]),
            'end': int(g[7]),
            'description': str(g[10]),
            'in_gene': bool(
                ((start <= int(g[6])) & (end > int(g[6]))) |
                ((start < int(g[7])) & (end >= int(g[7]))) |
                ((start >= int(g[6])) & (end <= int(g[7])))
            )
        } for g in genes
    ]
    return genes_arr


def collect_genes_for_regions(regions):
    """
    Collects the genes associated with each region. If there are multiple positions of interest within region,
    genes are found overlapping each position independently. If there are no positions of interest, genes are found
    overlapping the entire region. Found genes are attached to the region object and the region object is returned.
    """
    start = time.time()
    # for each region
    for region in regions:
        if 'top_positions' in region:
            region_genes = []
            for pos in region['top_positions']:
                genes_arr = get_genes(region['chromosome'], pos, pos)
                region_genes.append({
                    'i': int(pos),
                    'genes': genes_arr
                })
            region['genes'] = region_genes

        else:
            genes_arr = get_genes(
                region['chromosome'], 
                hilbert_pos_to_order(region['i'], region['order'], 14), 
                hilbert_pos_to_order(region['i'] + 1, region['order'], 14)
            )
            region['genes'] = genes_arr

    print(f"Time to collect gene information for regions: {time.time() - start}", flush=True)
    return regions


def collect_partial_paths_and_scores(
        regions,
        datatype='factor',
        factors_path=os.path.join(partial_csn_dir, 'paths/full_path_factors.max_order_resolution.h5'),
        scores_path=os.path.join(partial_csn_dir, 'scores/factor_scores.max_order_resolution.h5'),
        max_partial_order=13
    ):
    """
    Collects and returns partial path factors and scores for a list of regions. If region
    order > max partial path order, path factors/scores for max possible order are returned.
    """
    
    print(f"Collecting partial paths for {len(regions)} regions...", flush=True)
    start = time.time()

    file_path = factors_path if datatype == 'factor' else scores_path
    regions_data = []
    with h5py.File(file_path, 'r') as f:
        for r in regions:
            r_order = r['order']
            r_i = r['i']
            # check if region order is greater than max partial order
            if r_order > max_partial_order:
                r_i = hilbert_pos_to_order(r_i, r_order, max_partial_order)
                r_order = max_partial_order

            r_data = f[f"max_order_{r_order}.{r['chromosome']}"][r_i]
            regions_data.append(r_data)

    print(f"Time to collect partial paths for regions: {time.time() - start}", flush=True)

    return regions_data


def collect_full_paths_and_scores(
        regions,
        datatype='factor',
        factors_path=os.path.join(csn_dir, 'paths/full_path_factors.order_14_resolution.h5'),
        scores_path=os.path.join(csn_dir, 'scores/factor_scores.order_14_resolution.h5'),
    ):
    """
    Collects and returns full path factors and scores for a list of regions. If region
    order < 14, None is returned.
    """
    
    print(f"Collecting full paths for {len(regions)} regions...", flush=True)
    start = time.time()

    file_path = factors_path if datatype == 'factor' else scores_path
    regions_data = []
    with h5py.File(file_path, 'r') as f:
        for r in regions:
            r_order = r['order']
            r_i = r['i']

            if r_order == 14:
                r_data = f[f"{r['chromosome']}"][r_i]
                regions_data.append(r_data)
            else:
                regions_data.append(None)

    print(f"Time to collect full paths for regions: {time.time() - start}", flush=True)

    return regions_data


def generate_partial_path_return(regions, regions_factors, regions_scores):
    """
    Create return object for partial paths
    """
    for i in range(len(regions)):
        region = regions[i]
        region['path_factors'] = regions_factors[i].tolist()
        region['factor_scores'] = regions_scores[i].tolist()

    return json.dumps({"regions": regions})


@csns_bp.route("top_paths", methods=['POST'])
def top_paths():
    """
    Endpoint for fetching the top paths based on a set of filters and regions.
    """
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


@csns_bp.route("preview_filter", methods=['POST'])
def preview_filter():
    """
    Endpoint for previewing the number of regions a filter would provide at each order if selected
    """
    request_info = request.get_json()
    
    filters = request_info['filters']  # [{order, index, dataset_name}, ...]
    new_filter = request_info['newFilter']  # [{order, index}, ...]
    region = request_info['region']  # {order, chromosome, index}

    print(filters, new_filter)

    previews = get_factor_filter_previews(filters, new_filter, region)

    return_object = json.dumps({"preview_fractions": previews})
    return return_object


@csns_bp.route("top_paths_for_regions", methods=['POST'])
def native_order_resolution_regions():
    """
    Endpoint for fetching the top paths for a number of regions. The top N paths for each region are returned.
    """
    request_info = request.get_json()
    
    region_coordinates = request_info['regions']  # ['chr:start-end', ...] or [{chromosome, start, end}, ...]
    N = request_info['N'] if 'N' in request_info else 1 # number of paths to return per region (default 1)

    # convert region coordinates to regions object
    regions = [
        r if isinstance(r, dict) else {
            'chromosome': r.split(':')[0],
            'start': int(r.split(':')[1].split('-')[0]),
            'end': int(r.split(':')[1].split('-')[1])
        } for r in region_coordinates
    ]

    # get the scores for each region
    regions = get_scores_for_regions(regions, "path")

    # find the top N positions based on path score
    regions = get_top_scores_for_regions(regions, N)

    # collect factor scores for top positions
    regions = get_scores_for_regions(regions, "factor", by_position=True)

    # get the genes for each position in each region
    regions = collect_genes_for_regions(regions)

    # get the dehydrated paths for each position in each region
    regions = get_dehydrated_paths(regions=regions)

    # generate the return object
    return_object = generate_region_return(regions)
    
    return return_object


@csns_bp.route("paths_by_order", methods=['POST'])  # was partial_paths_for_regions
def paths_by_order():
    """
    Endpoint for collecting partial paths for a list of regions
    """
    request_info = request.get_json()
    
    regions = request_info['regions']  # [{chromosome, order, i}, ...]
    use_full_paths = request_info.get("use_full_paths", True)
    use_full_paths = str(use_full_paths).lower() == "true"

    # collect partial path factors and scores
    regions_partial_factors = collect_partial_paths_and_scores(regions, datatype='factor')
    regions_partial_scores = collect_partial_paths_and_scores(regions, datatype='score')
    
    # full and partial paths are stored in separate files
    if use_full_paths:
        # collect full path factors and scores for order 14 regions
        regions_full_factors = collect_full_paths_and_scores(regions, datatype='factor')
        regions_full_scores = collect_full_paths_and_scores(regions, datatype='score')

        # merge
        regions_factors = [
            regions_partial_factors[i] if r["order"] < 14 else regions_full_factors[i] for i, r in enumerate(regions)
        ]
        regions_scores = [
            regions_partial_scores[i] if r["order"] < 14 else regions_full_scores[i] for i, r in enumerate(regions)
        ]

    else:
        regions_factors = regions_partial_factors
        regions_scores = regions_partial_scores

    # get genes for each region
    regions = collect_genes_for_regions(regions)

    # generate the return object
    return_object = generate_partial_path_return(regions, regions_factors, regions_scores)

    return return_object