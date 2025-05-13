from flask import Blueprint, request
import numpy as np
import json
import h5py
import os
import time
from common import hilbert_pos_to_order, tabix_query
from constants import full_gene_path, partial_csn_dir, geneset_enrichment

# Create a Blueprint
csns_bp = Blueprint('csns_bp', __name__)


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


def generate_partial_path_return(regions, regions_factors, regions_scores, geneset_data):
    """
    Create return object for partial paths
    """
    for i in range(len(regions)):
        region = regions[i]
        region['path_factors'] = regions_factors[i].tolist()
        region['factor_scores'] = regions_scores[i].tolist()

    # add geneset object if present (not membership)
    return_object = {"regions": regions}
    if geneset_data is not None:
        return_object["genesets"] = geneset_data

    return json.dumps(return_object)


@csns_bp.route("paths_by_order", methods=['POST'])  # was partial_paths_for_regions
def paths_by_order():
    """
    Endpoint for collecting partial paths for a list of regions
    """
    request_info = request.get_json()
    
    regions = request_info['regions']  # [{chromosome, order, i}, ...]
    threshold = request_info.get("threshold", 0.05)  # default geneset p-value threshold is 0.05
    membership = str(request_info.get("membership", "False")).lower() == "true"  # default geneset membership is False, convert to boolean
    bonferroni = str(request_info.get("bonferroni", "False")).lower() == "true"  # default geneset bonferroni correction is False, convert to boolean

    # collect partial path factors and scores
    regions_factors = collect_partial_paths_and_scores(regions, datatype='factor')
    regions_scores = collect_partial_paths_and_scores(regions, datatype='score')

    # get genes for each region
    regions = collect_genes_for_regions(regions)

    if membership:
        # collect geneset membership for each region independently
        for region in regions:
            region_genes = [g["name"] for g in region["genes"]]
            p = geneset_enrichment.perform_geneset_hypergeometric_test(region_genes, membership=membership, bonferroni=bonferroni)
            geneset_data = geneset_enrichment.generate_geneset_return_data(p, threshold=threshold)
            region["genesets"] = geneset_data
    else:
        # collect geneset enrichments for region set as a whole
        genes_for_test = [g["name"] for region in regions for g in region["genes"]]
        p = geneset_enrichment.perform_geneset_hypergeometric_test(genes_for_test, membership=membership, bonferroni=bonferroni)
        geneset_data = geneset_enrichment.generate_geneset_return_data(p, threshold=threshold)

    # generate the return object
    return_object = generate_partial_path_return(regions, regions_factors, regions_scores, geneset_data if not membership else None)

    return return_object