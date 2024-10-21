from flask import Flask, request, abort, Response
from flask_cors import CORS
import numpy as np
from scipy.stats import hypergeom
import re
# import tabix
import sys
import os
from git import Repo
import json
# import h5py
from subprocess import Popen, PIPE
# from simsearch_functions import perform_simsearch, roiURL  # , load_per

# CSN
from csns import csns_bp
from data_filtering import data_filtering_bp
from gwas import gwas_bp
from geneset_enrichment import geneset_enrichment_bp


app = Flask(__name__)
CORS(app)
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024 # 16 MB


current_name = (__file__)
prod_name = '/data/narrativeGenomics/NarrativeGenomics/region_narration/region_narration_app/region-narration-flask/region_narration.py'

if current_name == prod_name:
    production = True
    server_root_dir = '/data/narrativeGenomics/NarrativeGenomics/region_narration'

    RN_VERSION = '0.0.1'
    RN_VERSION_STR = "{}/{}".format(Repo(server_root_dir, search_parent_directories=True).head.ref.name, RN_VERSION)

    # load search spaces for each scale
    data_dir = os.path.join(server_root_dir, 'data')

else:
    production = False
    RN_VERSION_STR = ''
    # load search spaces for each scale
    data_dir = '../../data/'

# print(production)

# Given a scale, build a path name for the search space
def build_search_space_path(scale, method='sliding_window'):
    if method == 'sliding_window':
        path = os.path.join(data_dir, 'sliding_window/regions_{}kb_CT20230427.DHS.factornormNMF.sum_enrich_corank.epilogos.statesum_rank.percentiles.limited_search_space_top_1000_no_overlap.npz'.format(scale))  # removed
    elif method == 'hilbert':
        path = os.path.join(data_dir, 'hilbert/hilbert_regions_order_{}_CT20230606.DHS.factornormNMF.sum_enrich_corank.epilogos.statesum_rank.percentiles.limited_search_space_top_1000_no_overlap.npz'.format(scale))  # removed
    elif method == 'hilbert_sfc':
        # path = os.path.join(data_dir, 'hilbert_sfc/hilbert_regions_order_{}_CT20230630.dhs.cs.sfc.ranks.limited_search_space_top_1000_no_overlap.npz'.format(scale))
        path = os.path.join(data_dir, 'hilbert_sfc/hilbert_regions_order_{}_CT20230718.dhs.cs.sfc.ranks.search_space.ft_2.rt_21.npz'.format(scale))  # removed
    return path


# Given a scale, build a path name for the search space
def build_geneset_bin_path(scale):
    path = os.path.join(data_dir, 'hilbert_sfc/hilbert_regions_order_{}.dhs.cs.sfc.ranks.search_space.ft_2.rt_21.mapped_extended_domain_GREATv4_gene_c5GOBP_selected_genesets.bin.compressed.CT20231023.npz'.format(scale))  # removed
    return path


# Given a scale, build a path name for the tabix query
def build_tabix_path(scale, method='sliding_window'):
    if method == 'sliding_window':
        filename = 'sliding_window/hg38.{}kb.regions.CT20230427.sorted.rownum4.factorranks5t38.bed.gz'.format(scale)  # removed
    elif method == 'hilbert':
        filename = 'hilbert/hilbert_order_{}_windows.sorted.rownum4.factorranks5t38.bed.gz'.format(scale)  # removed
    elif method == 'hilbert_sfc':
        filename = 'hilbert_sfc/20230718/hilbert_order_{}_windows.sorted.rownum4.sfc.ranks5t38.bed.gz'.format(scale)  # removed
    if production:
        path = 'http://explore.altius.org/tabix/ng/'
    else:
        path = 'https://resources.altius.org/~ctrader/public/genomic_narration/indexed/'
    return os.path.join(path, filename)


# Load search space data from a given scale
def load_ranks(search_space_path, method='sliding_window'):
    search_space = np.load(search_space_path)
    if ((method == 'sliding_window') | (method == 'hilbert')):
        ranks = search_space['percentiles']
    elif (method == 'hilbert_sfc'):
        ranks = search_space['ranks']
    coords = search_space['coordinates']
    return ranks, coords


# Load search space data from a given scale
def load_geneset_bin(geneset_bin_path):
    geneset_bin = np.load(geneset_bin_path)
    search_space_geneset_bin = geneset_bin['search_space_geneset_bin']
    geneset_order = geneset_bin['geneset_order']
    return search_space_geneset_bin, geneset_order


# load search spaces for all scales into dict
# sliding_window_seach_spaces = {}
# hilbert_search_spaces = {}
hilbert_sfc_search_spaces = {}
hilbert_sfc_geneset_bin = {}
# for scale in [1, 2, 5, 10, 15, 25, 50, 75, 100]:
#     region_method = 'sliding_window'
#     sliding_window_seach_spaces[scale] = load_ranks(build_search_space_path(scale, region_method), region_method)
# for scale in [4, 5, 6, 7, 8, 9, 10, 11]:
#     region_method = 'hilbert'
#     hilbert_search_spaces[scale] = load_ranks(build_search_space_path(scale, region_method), region_method)

# for scale in [4, 5, 6, 7, 8, 9, 10, 11]:
#     region_method = 'hilbert_sfc'
#     hilbert_sfc_search_spaces[scale] = load_ranks(build_search_space_path(scale, region_method), region_method)
#     hilbert_sfc_geneset_bin[scale] = load_geneset_bin(build_geneset_bin_path(scale))

# Functions for narration
# run tabix query
def tabix_query(filename, query):
    # Use subprocess with tabix and return hits
    process = Popen(['tabix', '-f', filename, query], stdout=PIPE)

    hits = []
    for line in process.stdout:
        hits.append(np.char.decode(line.strip().split()))
    hits = np.array(hits)
    return hits


# determine roi from query coordinates
def load_roi(genomic_location, link_path):
    # determine chromosome, start, and end from query
    query_start_end = genomic_location.split(':')[1].split('-')
    query_start = int(query_start_end[0])
    query_end = int(query_start_end[1])

    # retrieve possible rois
    possible_coords = tabix_query(link_path, genomic_location)

    # find closest roi
    start_dist = np.abs(possible_coords[:, 1].astype(int) - query_start)
    end_dist = np.abs(possible_coords[:, 2].astype(int) - query_end)
    total_dist = start_dist + end_dist
    min_dist = np.min(total_dist)
    min_dist_region_inds = np.where(total_dist == min_dist)[0]
    middle_ind = min_dist_region_inds[np.floor(len(min_dist_region_inds) / 2).astype(int)]
    query_region = possible_coords[middle_ind]

    roi = int(query_region[3])
    roi_ranks = query_region[4:].astype('float32')
    roi_coords = query_region[:3]

    return roi, roi_ranks, roi_coords


# choose search space from scale options
def choose_search_space(scale, method='sliding_window'):
    if method == 'sliding_window':
        ranks, coords = sliding_window_seach_spaces[scale]
    elif method == 'hilbert':
        ranks, coords = hilbert_search_spaces[scale]
    elif method == 'hilbert_sfc':
        ranks, coords = hilbert_sfc_search_spaces[scale]
    return ranks, coords

# limit ranks to just what the user wants to include
def filter_factors_in_ranks(included_factors_str, ranks, roi_ranks):
    metric_factors = np.int_(
        included_factors_str.strip('[').strip(']').split(',')
    )
    ranks = ranks[:, metric_factors]
    roi_ranks = roi_ranks[metric_factors]

    return ranks, roi_ranks

# find the factors we want to search over
def refine_factors_to_search_over(query_factor_thresh, roi_ranks):
    query_factors_above_thresh = np.where(roi_ranks > query_factor_thresh)[0]
    # if no factors above threshold, search best one
    if len(query_factors_above_thresh) == 0:
        sort_order = np.argsort(-roi_ranks)
        included_factors = sort_order[:1]
    else:
        included_factors = query_factors_above_thresh
    return included_factors

# # given our roi, find the similar regions
# def find_most_similar_region_indices(roi_ranks, ranks, roi, roi_coords, coords, included_factors, num_similar_regions):
#     (top_regions, top_region_coords, top_regions_sim,
#      top_regions_per, query_factor_order,
#      target_factors, pval,
#      Z, query_sim_prod_ranks) = perform_simsearch(roi_ranks, ranks, roi, roi_coords, coords, included_factors, num_similar_regions)
    
#     return top_regions, top_region_coords, top_regions_sim, top_regions_per, query_factor_order, target_factors, pval, Z, query_sim_prod_ranks


# # generate roiURL to view similar regions in index browser
# def setup_roiURLs(get_roiURL, target_factors, included_factors, roi, top_region_coords):
#     roiURLs = []
#     for dl in range(1, len(included_factors) + 1):
#         if get_roiURL == 'Tru':  # 'True'
#             # roiURLs.append('roiURL not available')
#             roiURLs.append(
#                 roiURL(
#                     roi, 
#                     top_region_coords[dl-1],
#                     dl, 
#                     '_'.join(target_factors.astype(str))
#                 )
#             )
#         else:
#             roiURLs.append('roiURL not available')
#     return roiURLs


# # reorder data based on ranks in query region/roi
# def reorder_by_query_region(top_regions_per, query_factor_order, pval, Z):
#     return_top_regions_per = top_regions_per[:, :, query_factor_order]
#     return_pval = pval[:, query_factor_order]
#     return_Z = Z[:, query_factor_order]
#     return return_top_regions_per, return_pval, return_Z


# # build array of objects for roi and each similar region at each detail level
# def generate_simsearch_result_array(num_similar_regions, return_top_regions_per, included_factors, top_regions, top_regions_sim, top_region_coords, query_sim_prod_ranks, query_factor_order, target_factors, return_Z, return_pval, roiURLs):
#     return_arr_full = []
#     # if the best factor for the query region has a non-zero rank
#     if return_top_regions_per[0, 0, 0] > 0:
#         for dl in range(len(included_factors)):
#             # for the roi, return results as json array
#             # each element is an object of scores for each metric
#             return_arr = []
#             for rank in range(num_similar_regions + 1):
#                 # if similar region has non-zero score
#                 if float(top_regions_sim[dl, rank]) > 0:
#                     coord = top_region_coords[dl, rank]
#                     region_dict = {
#                         'rank': rank,
#                         'coordinates': coord[0] + ':' + coord[1] + '-' + coord[2],
#                         'region_ind': int(top_regions[dl, rank]),
#                         'score': float(top_regions_sim[dl, rank]),
#                         'percentiles': return_top_regions_per[dl, rank].tolist(),
#                         'prod_rank': query_sim_prod_ranks[dl, rank]
#                     }
#                     if rank == 0:
#                         region_dict['region_ind'] = 'NA'
#                         region_dict['score'] = 'NA'
#                         region_dict['factor_order'] = query_factor_order.tolist()
#                         region_dict['target_factors'] = target_factors.tolist()
#                         region_dict['top_region_factor_Z_score'] = return_Z[dl].tolist()
#                         region_dict['top_region_factor_pval'] = return_pval[dl].tolist()
#                         region_dict['roiURL'] = roiURLs[dl]
#                         region_dict['detail_level'] = dl + 1
#                         # region_dict['query_prod_rank'] = query_sim_prod_ranks[dl]
#                     return_arr.append(region_dict)
#             return_arr_full.append(return_arr)
#     return return_arr_full

# build array of objects for roi and each similar region at each detail level
def generate_simsearch_result_array(num_similar_regions, top_regions_ranks, selected_factors, top_regions, top_regions_score, top_regions_coords):
    # if the best factor for the query region has a non-zero rank
    if np.sum(top_regions_ranks[0]) > 0:
        # for the roi, return results as json array
        # each element is an object of scores for each metric
        return_arr = []
        for rank in range(num_similar_regions + 1):
            # if similar region has non-zero score
            if float(top_regions_score[rank]) > 0:
                coord = top_regions_coords[rank]
                region_dict = {
                    'rank': rank,
                    'coordinates': coord[0] + ':' + coord[1] + '-' + coord[2],
                    'region_ind': int(top_regions[rank]),
                    'score': float(top_regions_score[rank]),
                    'percentiles': top_regions_ranks[rank].tolist(),
                    # 'prod_rank': query_sim_prod_ranks[dl, rank]
                }
                if rank == 0:
                    region_dict['region_ind'] = 'NA'
                    region_dict['score'] = 'NA'
                    region_dict['selected_factors'] = selected_factors.tolist()
                    # region_dict['factor_order'] = query_factor_order.tolist()
                    # region_dict['target_factors'] = target_factors.tolist()
                    # region_dict['top_region_factor_Z_score'] = return_Z[dl].tolist()
                    # region_dict['top_region_factor_pval'] = return_pval[dl].tolist()
                    # region_dict['roiURL'] = roiURLs[dl]
                    # region_dict['detail_level'] = dl + 1
                    # region_dict['query_prod_rank'] = query_sim_prod_ranks[dl]
                return_arr.append(region_dict)
    return return_arr


# find the factors we want to search over
def setup_selected_factors(selected_factors_str='', roi_ranks=[], query_factor_thresh=0.5):
    if selected_factors_str != '':
        selected_factors = np.int_(
            selected_factors_str.strip('[').strip(']').split(',')
        )
    elif len(roi_ranks) > 0:
        selected_factors = refine_factors_to_search_over(query_factor_thresh, roi_ranks)
    return selected_factors


def determine_factor_order(top_regions, selected_factors, ranks):
    top_region_ranks = ranks[top_regions]

    # order selected factors first
    selected_factors_best_order = selected_factors[
        np.argsort(
            -np.sum(
                top_region_ranks[:, selected_factors],
                axis=0
            )
        )
    ]

    # order remaining factors
    regions_to_consider = np.where(
        np.count_nonzero(
            ranks[:, selected_factors_best_order],
            axis=1
        ) > 0
    )[0]

    total_factors = np.arange(ranks.shape[1])
    leftover_factors = total_factors[~np.isin(total_factors, selected_factors_best_order)]
    leftover_factors_best_order = leftover_factors[
        np.argsort(
            -np.sum(
                ranks[regions_to_consider][:, leftover_factors],
                axis=0
            ),
            kind='stable'
        )
    ]

    factor_order = np.append(selected_factors_best_order, leftover_factors_best_order)
    return factor_order


# search for the best regions for selected factors
def search_over_selected_factors(num_regions, ranks, coords, selected_factors, convolution=[], roi_coords=[], weight=0.8):
    # fusion
    if len(convolution) > 0:
        data_for_fusion = convolution[:, selected_factors]
    else:
        data_for_fusion = ranks[:, selected_factors]
    selected_factor_search_space_score = 1 - np.prod(1 - weight * data_for_fusion, axis=1)

    # get the best regions
    if len(roi_coords) > 0:
        # filter out regions that overlap with the query region
        mask = np.where(~
            (
                (coords[:, 0] == roi_coords[0]) & 
                (
                    (
                        (coords[:, 1].astype(int) > int(roi_coords[1])) &
                        (coords[:, 1].astype(int) < int(roi_coords[2]))
                    ) |
                    (
                        (coords[:, 2].astype(int) > int(roi_coords[1])) &
                        (coords[:, 2].astype(int) < int(roi_coords[2]))
                    ) |
                    (
                        (coords[:, 1].astype(int) == int(roi_coords[1])) &
                        (coords[:, 2].astype(int) == int(roi_coords[2]))
                    )
                )
            )
        )[0]
        top_regions = mask[np.argsort(-selected_factor_search_space_score[mask])][:num_regions]
    else:
        top_regions = np.argsort(-selected_factor_search_space_score)[:num_regions]

    # factor order
    # factor_order = determine_factor_order(top_regions, selected_factors, ranks)

    # find info for each top region
    top_regions_ranks = ranks[top_regions]#[:, factor_order]
    top_regions_coords = coords[top_regions]
    top_regions_score = selected_factor_search_space_score[top_regions]
    return top_regions, top_regions_ranks, top_regions_coords, top_regions_score#, factor_order


# build array of objects for each region at each detail level
def generate_search_by_factor_result_array(top_regions, top_regions_coords, top_regions_score, top_regions_ranks, selected_factors): # factor_order
    # for each region, return relevant information
    search_by_factor_return_arr = []
    for rank in range(len(top_regions)):
        coord = top_regions_coords[rank]
        # check for non-zero score
        if float(top_regions_score[rank]) > 0:
            region_dict = {
                'rank': rank,
                'coordinates': coord[0] + ':' + coord[1] + '-' + coord[2],
                'region_ind': int(top_regions[rank]),
                'score': float(top_regions_score[rank]),
                'percentiles': top_regions_ranks[rank].tolist()
            }
            if rank == 0:
                # region_dict['factor_order'] = factor_order.tolist()
                region_dict['selected_factors'] = selected_factors.tolist()
            search_by_factor_return_arr.append(region_dict)
    return search_by_factor_return_arr

def generate_narration_result_array(roi, roi_ranks, roi_coords):
    narration_return_arr = [{'roi': roi, 'percentiles': roi_ranks.tolist(), 'coordinates': roi_coords.tolist()}]
    return narration_return_arr

def calculate_convolution(roi_ranks, ranks):
    conv = ranks * roi_ranks
    return conv

def append_roi_information(top_regions, top_regions_ranks, top_regions_coords, top_regions_score, roi, roi_ranks, roi_coords): #factor_order
    appended_top_regions = np.append(roi, top_regions)

    appended_top_regions_ranks = np.zeros((top_regions_ranks.shape[0] + 1, top_regions_ranks.shape[1]))
    appended_top_regions_ranks[0]  = roi_ranks#[factor_order]
    appended_top_regions_ranks[1:] = top_regions_ranks

    appended_top_regions_coords = np.zeros((top_regions_coords.shape[0] + 1, top_regions_coords.shape[1]), dtype=object)
    appended_top_regions_coords[0]  = roi_coords
    appended_top_regions_coords[1:] = top_regions_coords

    appended_top_regions_score = np.append(1, top_regions_score)

    return appended_top_regions, appended_top_regions_ranks, appended_top_regions_coords, appended_top_regions_score


def geneset_hypergeom(region_inds, N, M, n, geneset_order, search_space_geneset_bin):
    k = np.sum(search_space_geneset_bin[region_inds], axis=0)
    inds_to_test = np.where(k)[0]

    test = hypergeom.sf(k[inds_to_test], M, n[inds_to_test], N)
    tested_genesets = geneset_order[inds_to_test]
    test = test * len(geneset_order)  # bonferroni correction
    test = np.where(test > 1, 1, test)
    return tested_genesets, inds_to_test, test


def parse_geneset_enrichment_query_data(query_data):
    query_sim_region_inds = []
    for i in range(len(query_data)):
        query_sim_region_inds.append(query_data[i]['region_ind'])
    
    return query_sim_region_inds


@app.route('/', methods=['GET'])
def about():
    python_version = sys.version
    return "<html><body><div>region-narration ({})</div><div>Python {}</div></body></html>".format(RN_VERSION_STR, python_version)


# find the top regions for selected factors
@app.route('/search_by_factor', methods=['POST'])
def search_by_factor():
    request_info = request.get_json()
    scale = request_info['scale']
    selected_factors_str = request_info['factors']
    num_regions = int(request_info['numRegions'])
    region_method = request_info['regionMethod']
    weight = float(request_info['fusionWeight'])
    print((scale, selected_factors_str, num_regions, region_method, weight), flush=True)

    # choose search space
    ranks, coords = choose_search_space(scale, region_method)

    # get the selected factors
    selected_factors = setup_selected_factors(selected_factors_str)

    # get the top regions for the selected factors
    top_regions, top_regions_ranks, top_regions_coords, top_regions_score = search_over_selected_factors(num_regions, ranks, coords, selected_factors, weight=weight)

    # build return array
    search_by_factor_return_arr = generate_search_by_factor_result_array(top_regions, top_regions_coords, top_regions_score, top_regions_ranks, selected_factors)

    # return results
    json_return = json.dumps(search_by_factor_return_arr)
    return json_return


@app.route('/simsearch', methods=['POST'])
def simsearch():
    request_info = request.get_json()
    genomic_location = request_info['location']
    scale = request_info['scale']
    included_factors_str = request_info['includedFactors']
    selected_factors_str = request_info['selectedFactors']
    get_roiURL = request_info['roiURL']
    query_factor_thresh = float(request_info['queryFactorThresh'])
    num_similar_regions = int(request_info['numSimilarRegions'])
    region_method = request_info['regionMethod']
    weight = float(request_info['fusionWeight'])

    # get roi
    roi, roi_ranks, roi_coords = load_roi(genomic_location, build_tabix_path(scale, region_method))
    print(build_tabix_path(scale, region_method), flush=True)

    # choose search space
    ranks, coords = choose_search_space(scale, region_method)

    # limit ranks to just the factors included with included_factors_str
    ranks, roi_ranks = filter_factors_in_ranks(included_factors_str, ranks, roi_ranks)
    if len(roi_ranks[roi_ranks != 0]) == 0:
        abort(404)

    # get the selected factors
    selected_factors = setup_selected_factors(selected_factors_str, roi_ranks, query_factor_thresh)

    # calculate convolution scores between roi and search space
    convolution = calculate_convolution(roi_ranks, ranks)

    top_regions, top_regions_ranks, top_regions_coords, top_regions_score = search_over_selected_factors(num_similar_regions, ranks, coords, selected_factors, convolution, roi_coords=roi_coords, weight=weight)

    top_regions, top_regions_ranks, top_regions_coords, top_regions_score = append_roi_information(top_regions, top_regions_ranks, top_regions_coords, top_regions_score, roi, roi_ranks, roi_coords)

    # generate result array
    return_arr = generate_simsearch_result_array(num_similar_regions, top_regions_ranks, selected_factors, top_regions, top_regions_score, top_regions_coords) #query_sim_prod_ranks

    return return_arr




# def simsearch_old():
#     request_info = request.get_json()
#     genomic_location = request_info['location']
#     scale = request_info['scale']
#     included_factors_str = request_info['factors']
#     get_roiURL = request_info['roiURL']
#     query_factor_thresh = float(request_info['queryFactorThresh'])
#     num_similar_regions = int(request_info['numSimilarRegions'])
#     region_method = request_info['regionMethod']

#     # get roi
#     roi, roi_ranks, roi_coords = load_roi(genomic_location, build_tabix_path(scale, region_method))
#     print(build_tabix_path(scale, region_method), flush=True)

#     # choose search space
#     ranks, coords = choose_search_space(scale, region_method)

#     # limit ranks to just the factors included with included_factors_str
#     ranks, roi_ranks = filter_factors_in_ranks(included_factors_str, ranks, roi_ranks)

#     # find the factors we actually want to search over/the user would care about
#     included_factors = refine_factors_to_search_over(query_factor_thresh, roi_ranks)

#     # find most similar region inds
#     top_regions, top_region_coords, top_regions_sim, top_regions_per, query_factor_order, target_factors, pval, Z, query_sim_prod_ranks = find_most_similar_region_indices(roi_ranks, ranks, roi, roi_coords, coords, included_factors, num_similar_regions)

#     # dual track roiURL
#     roiURLs = setup_roiURLs(get_roiURL, target_factors, included_factors, roi, top_region_coords)

#     # reorder based on query region
#     return_top_regions_per, return_pval, return_Z = reorder_by_query_region(top_regions_per, query_factor_order, pval, Z)

#     # generate result array
#     return_arr_full = generate_simsearch_result_array(num_similar_regions, return_top_regions_per, included_factors, top_regions, top_regions_sim, top_region_coords, query_sim_prod_ranks, query_factor_order, target_factors, return_Z, return_pval, roiURLs)

#     # result array as json
#     json_return = json.dumps(return_arr_full)

#     return json_return


# find the top regions for selected factors
@app.route('/narration', methods=['POST'])
def narration():
    request_info = request.get_json()
    genomic_location = request_info['location']
    scale = request_info['scale']
    region_method = request_info['regionMethod']

    # get roi
    roi, roi_ranks, roi_coords = load_roi(genomic_location, build_tabix_path(scale, region_method))

    # generate result array
    narration_return_arr = generate_narration_result_array(roi, roi_ranks, roi_coords)

    # result array as json
    json_return = json.dumps(narration_return_arr)

    return json_return


# perform a geneset enrichment test for desired regions
@app.route('/geneset_enrichment', methods=['POST'])
def geneset_enrichment():
    request_info = request.get_json()
    query_data = request_info['data']
    scale = request_info['scale']


    # get the search space x geneset matrix for the correct order
    search_space_geneset_bin, geneset_order = hilbert_sfc_geneset_bin[scale]

    # parse the data to get the region indices and the region coordinates
    query_sim_region_inds = parse_geneset_enrichment_query_data(query_data)
        
    # perform hypergeometric test for enrichment
    num_sim_regions = len(query_sim_region_inds)
    num_regions = search_space_geneset_bin.shape[0]
    gs_region_counts = np.sum(search_space_geneset_bin, axis=0)

    query_tested_gs, query_tested_gs_inds, query_tested_gs_pval = geneset_hypergeom(
        region_inds=query_sim_region_inds, 
        N=num_sim_regions,
        M=num_regions, 
        n=gs_region_counts, 
        geneset_order=geneset_order,
        search_space_geneset_bin=search_space_geneset_bin
    )

    return_geneset_enrichment = []
    for gs, i, p in zip(query_tested_gs, query_tested_gs_inds, query_tested_gs_pval):
        return_geneset_enrichment.append({
            'geneset_name': gs, 
            'geneset_index': int(i), 
            'geneset_pval': float(p)
        })
        
    json_return = json.dumps(return_geneset_enrichment)

    return json_return


# =======================================================================================
# YC240401 Shared factor anaylsis for UMAP
# =======================================================================================

import numpy as np
import pandas as pd
import h5py

# GilbertSFCH5_DIR = '/home/ychoi/projects/hilbert-genome/data/20240328_GilbertDataInH5'
GilbertSFCH5_DIR = '/data/narrativeGenomics/NarrativeGenomics/region_narration/data/hilbert_sfc_h5/20240328_GilbertDataInH5'

# Map chromosome location to index at specific Gilbert order
def map_to_order_indice(chr_str, order, df_indices):
    # Parse chr_str
    chr,start_end = chr_str.split(':')
    start,end = [int(i) for i in start_end.split('-')]
    # Get index 
    segment_size = 4**(14-order)
    chr_index = start // segment_size
    chr_offset = df_indices[str(order)][chr] 
    index = chr_index + chr_offset
    return index

# Use fancy indexing on h5
def fancy_indexing(hf, index):
    # Create unique ordered version of index for hf
    index_sorted = np.sort(np.unique(index)) # index for hf
    hfidx2idx = dict(zip(index_sorted,np.arange(len(index_sorted))))
    # Get array from hf then reindex
    arr_hf = np.array(hf['arr'][index_sorted],dtype='float32')
    arr = arr_hf[[hfidx2idx[i] for i in index]]
    return arr

# Given dataset and order, create path to h5 file
def get_h5_path(ds,order):
    basedir = GilbertSFCH5_DIR
    path = f"{basedir}/{ds}/offset_0/order_{order}_windows_annotation_sfc_scores.h5"
    return path

# Given variable metadata for factors and kwargs for columns in the metadata, return boolean index
def subset_by_var(var, **kwargs):
    var_idx = np.repeat(True, len(var))
    for k,v in kwargs.items():
        var_idx = np.all([var_idx,(var[k] == v)],axis=0)
    return var_idx

# Given SFC array of (sample,feature), return top n arg per sample (sample, n)
def argnmax(arr, n=1, axis=1, ordered=True):
    if ordered:
        idx = np.argsort(-arr, axis=axis)[:,:n]
    else:
        idx = np.argpartition(arr, -n, axis=axis)[:,-n:]
    return idx

# Given SFC array of (sample,feature), return top args per sample and index of features considered
def get_top_factor(arr, var, n=1, ordered=True, exclude_LAD=False, **kwargs):
    """
    """
    var_idx = subset_by_var(var, **kwargs)
    if exclude_LAD:
        var_idx[subset_by_var(var,dataset='LAD')] = False
    arr_keep = arr[:,var_idx]
    topargs = argnmax(arr_keep,n=10,ordered=ordered)
    zero_idx = np.array([arr_keep[i,topargs[i]] for i in range(len(topargs))]) == 0
    topargs[zero_idx] = -1
    return topargs, var_idx

def shared_factor_percentage(top_factors, n_shared=10):
    """
    per factor, show percentage of the samples the factor was top factor in
    """
    top_shared_factor_counts = pd.Series(top_factors.reshape(-1)).value_counts()
    top_shared_factor_counts_per_sample = top_shared_factor_counts / len(top_factors)
    return top_shared_factor_counts_per_sample

def perform_shared_factor_analysis(arr, var, n=10, min_percentage=0.2):
    """
    return
        (mean shared factor percentage, 
         list of shared factor names in decreasing order,
         individual shared factor perecentage)
    """
    top_factors,var_idx = get_top_factor(arr,n=n,var=var)
    shared_factor = shared_factor_percentage(top_factors)
    shared_factor = shared_factor[shared_factor.index != -1]
    score = shared_factor[:n].sum() / n
    shared_factor_name = var['scaled factor'][var_idx].values[shared_factor.index]
    shared_factor.index = shared_factor_name
    shared_factor = shared_factor[shared_factor.values > min_percentage]
    gil = [score,shared_factor_name,shared_factor]
    return gil


class ScaledGilbertDataLoader():
    """
    class for loading cross scale Gilbert SFC data from h5 files
    """
#
    def __init__(self):
#
        # Preload h5 files
        self.datasets = ['DHS','chromatin_states','TF']
        self.orders = np.arange(4,13)
#
        h5_dict = {}
        for ds in self.datasets:
            for order in self.orders:
                h5path = get_h5_path(ds,order)
                h5_dict[f"{ds}-{order}"] = h5py.File(h5path,'r')
#
        self.h5_dict = h5_dict
#
        # Load indices where each chromosome begins
        FILE_DIR = f'{GilbertSFCH5_DIR}/chromosome_indices.csv'
        self.df_indices = pd.read_csv(FILE_DIR,index_col=0)
#
        # Load var, variable metadata
        self._generate_var()
#
    def load_segment(self, chr,start,end):
        arrs = []
        for ds in self.datasets: 
            for order in self.orders:
                idx = map_to_order_indice(chr,start,end,order)
                hf = self.h5_dict[f"{ds}-{order}"]
                arr_ = hf['arr'][idx]
                arrs.append(arr_)
        arrs = np.concatenate(arrs,axis=1)
        return arrs
#
    def load_segments(self, chr_strs):
        arrs = []
        for ds in self.datasets: 
            for order in self.orders:
                idx = [map_to_order_indice(chr_str=chr_str,order=order,df_indices=self.df_indices) for chr_str in chr_strs]
                hf = self.h5_dict[f"{ds}-{order}"]
                arr_ = fancy_indexing(hf,idx)
                arrs.append(arr_)
        arrs = np.concatenate(arrs,axis=1)
        return arrs
#
    def _generate_var(self):
        var = []
        for sds,h5 in self.h5_dict.items():
            ds,order = sds.split('-')
            factors = h5['columns'][:].astype('str')
            df_var = pd.DataFrame({'Scaled Dataset': sds,
                                    'dataset': ds,
                                    'order': order, 
                                    'factor': factors,
                                    'scaled factor': pd.Series(factors) + '-' + order,
                                    })
            var.append(df_var)
        self.var = pd.concat(var)

class ComponentNarrationLoader():
    #
    def __init__(self, model_dir):
        self.model_dir = model_dir
        self.load_component(model_dir)
    #
        # Load indices where each chromosome begins
        FILE_DIR = f'{model_dir}/chromosome_indices.csv'
        self.df_indices = pd.read_csv(FILE_DIR,index_col=0)
    #
    def load_component(self,model_dir):
        """
        Input
            model_dir: directory containing the model with whole genome components and component narration
        """
        h5dir = os.path.join(model_dir,'component.h5')
        comp_meta_dir = os.path.join(model_dir,'CompAnnotFull.csv')
        dom_dir = os.path.join(model_dir,f'top_10_component_per_segment-raw.h5')
        comp_narr_full = pd.read_csv(comp_meta_dir,index_col=0)
        self.comp_narr_full = comp_narr_full
        hf = h5py.File(h5dir,'r')
        self.comp_hf = hf['components']
        hf = h5py.File(dom_dir,'r')
        self.dom_hf = hf['component_index']
    #
    def load_segments(self, chr_strs, order=12):
        idx = [map_to_order_indice(chr_str=chr_str,order=order,df_indices=self.df_indices) for chr_str in chr_strs]
        hf = self.dom_hf
        dom = fancy_indexing(hf,idx)
        return dom
    #    
    def narrate_region(self, chr_strs, min_support=0.2, max_n_comp=5):
        dom = self.load_segments(chr_strs)
        narr = pd.Series(dom[:,:5].flatten()).value_counts()
        narr = narr.reset_index()
        narr.columns = ['CompIndex','Count']
        narr['CompNarr'] = self.comp_narr_full.FullNarration.values[narr.CompIndex]
        narr['Support']  = narr.Count / len(chr_strs)
        narr = narr.iloc[:max_n_comp]
        narr = narr[narr.Support >= min_support]
        return narr


sgdl = ScaledGilbertDataLoader()
cnl  = ComponentNarrationLoader(GilbertSFCH5_DIR)

# find the top regions for selected factors
@app.route('/get_shared_factor', methods=['POST'])
def get_shared_factor():
    request_info = request.get_json()
    chr_strs = request_info['chr_strs']
    arr = sgdl.load_segments(chr_strs)
    results = perform_shared_factor_analysis(arr=arr,var=sgdl.var)
    json_return = json.dumps({
        "SharedFactor": list(results[2].index.values),
        "SharedFactorPercentage": list(results[2].values),
        "MeanSharedFactorPercentage": results[2].mean(),
        }
    )
    return json_return

@app.route('/narrate_regions', methods=['POST'])
def narrate_regions():
    request_info = request.get_json()
    chr_strs = request_info['chr_strs']
    df_narr = cnl(chr_strs)
    json_return = json.dumps({
        "CopmonentNarration": list(df_narr.CompNarr.values),
        "Support": list(df_narr.Support.values),
        }
    )
    return json_return





##########CSN#########
app.register_blueprint(csns_bp, url_prefix='/api/csns') # Register the Blueprint at /api/csns
#######Filtering######
app.register_blueprint(data_filtering_bp, url_prefix='/api/dataFiltering') # Register the Blueprint at /api/dataFiltering
#########GWAS#########
app.register_blueprint(gwas_bp, url_prefix='/api/gwas') # Register the Blueprint at /api/gwas
##Geneset Enrichment##
app.register_blueprint(geneset_enrichment_bp, url_prefix='/api/genesetEnrichment') # Register the Blueprint at /api/genesetEnrichment




if __name__ == '__main__':
    app.run(host='0.0.0.0')
