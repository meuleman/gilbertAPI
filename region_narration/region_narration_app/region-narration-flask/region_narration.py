from flask import Flask, request
from flask_cors import CORS
import numpy as np
# import tabix
import sys
import os
from git import Repo
import json
# import h5py
from subprocess import Popen, PIPE
from simsearch_functions import simsearch, roiURL  # , load_per

app = Flask(__name__)
CORS(app)


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

print(production)

# Given a scale, build a path name for the search space
def build_search_space_path(scale, method='sliding_window'):
    if method == 'sliding_window':
        path = os.path.join(data_dir, 'sliding_window/regions_{}kb_CT20230427.DHS.factornormNMF.sum_enrich_corank.epilogos.statesum_rank.percentiles.limited_search_space_top_1000_no_overlap.npz'.format(scale))
    elif method == 'hilbert':
        path = os.path.join(data_dir, 'hilbert/hilbert_regions_order_{}_CT20230606.DHS.factornormNMF.sum_enrich_corank.epilogos.statesum_rank.percentiles.limited_search_space_top_1000_no_overlap.npz'.format(scale))
    elif method == 'hilbert_sfc':
        # path = os.path.join(data_dir, 'hilbert_sfc/hilbert_regions_order_{}_CT20230630.dhs.cs.sfc.ranks.limited_search_space_top_1000_no_overlap.npz'.format(scale))
        path = os.path.join(data_dir, 'hilbert_sfc/hilbert_regions_order_{}_CT20230718.dhs.cs.sfc.ranks.search_space.ft_2.rt_21.npz'.format(scale))
    return path


# Given a scale, build a path name for the tabix query
def build_tabix_path(scale, method='sliding_window'):
    if method == 'sliding_window':
        filename = 'sliding_window/hg38.{}kb.regions.CT20230427.sorted.rownum4.factorranks5t38.bed.gz'.format(scale)
    elif method == 'hilbert':
        filename = 'hilbert/hilbert_order_{}_windows.sorted.rownum4.factorranks5t38.bed.gz'.format(scale)
    elif method == 'hilbert_sfc':
        filename = 'hilbert_sfc/20230718/hilbert_order_{}_windows.sorted.rownum4.sfc.ranks5t38.bed.gz'.format(scale)
    if production:
        path = 'http://explore.altius.org/tabix/ng/'
    else:
        path = 'https://resources.altius.org/~ctrader/public/genomic_narration/indexed/'
    return os.path.join(path, filename)


# Load search space data from a given scale
def load_ranks(search_space_path, method='sliding_window'):
    search_space = np.load(search_space_path)
    if ((method == 'sliding_window') | (method == 'hilbert')):
        percents = search_space['percentiles']
    elif (method == 'hilbert_sfc'):
        percents = search_space['ranks']
    coords = search_space['coordinates']
    return percents, coords


# load search spaces for all scales into dict
sliding_window_seach_spaces = {}
hilbert_search_spaces = {}
hilbert_sfc_search_spaces = {}
for scale in [1, 2, 5, 10, 15, 25, 50, 75, 100]:
    region_method = 'sliding_window'
    sliding_window_seach_spaces[scale] = load_ranks(build_search_space_path(scale, region_method), region_method)
for scale in [4, 5, 6, 7, 8, 9, 10, 11]:
    region_method = 'hilbert'
    hilbert_search_spaces[scale] = load_ranks(build_search_space_path(scale, region_method), region_method)
for scale in [4, 5, 6, 7, 8, 9, 10, 11]:
    region_method = 'hilbert_sfc'
    hilbert_sfc_search_spaces[scale] = load_ranks(build_search_space_path(scale, region_method), region_method)


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
    global roi, roi_percents, roi_coords

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
    roi_percents = query_region[4:].astype('float32')
    roi_coords = query_region[:3]

    return


# choose search space from scale options
def choose_search_space(scale, method='sliding_window'):
    global percents, coords
    if method == 'sliding_window':
        percents, coords = sliding_window_seach_spaces[scale]
    elif method == 'hilbert':
        percents, coords = hilbert_search_spaces[scale]
    elif method == 'hilbert_sfc':
        percents, coords = hilbert_sfc_search_spaces[scale]

# find the factors we want to search over
def setup_included_factors(included_factors_str, query_factor_thresh):
    global included_factors, metric_factors, percents, roi_percents

    metric_factors = np.int_(
        included_factors_str.strip('[').strip(']').split(',')
    )
    # limit ranks ro just metric_factors
    # print(percents, flush=True)
    percents = percents[:, metric_factors]
    roi_percents = roi_percents[metric_factors]

    query_factors_above_thresh = np.where(roi_percents > query_factor_thresh)[0]

    # if no factors above threshold, search best one
    if len(query_factors_above_thresh) == 0:
        sort_order = np.argsort(-roi_percents)
        # included_factors = sort_order[
        #     np.isin(sort_order, metric_factors)
        # ][:1]
        included_factors = sort_order[:1]
    else:
        # included_factors = metric_factors[
        #     np.isin(metric_factors, query_factors_above_thresh)
        # ]
        included_factors = query_factors_above_thresh
    # included_factors=np.arange(34)

# given our roi, find the similar regions
def find_most_similar_region_indices(num_similar_regions):
    global top_regions, top_region_coords, top_regions_sim, top_regions_per, query_factor_order, target_factors, pval, Z, query_sim_prod_ranks
    (top_regions, top_region_coords, top_regions_sim,
     top_regions_per, query_factor_order,
     target_factors, pval,
     Z, query_sim_prod_ranks) = simsearch(roi_percents, percents, roi, roi_coords, coords, included_factors, num_similar_regions)


# generate roiURL to view similar regions in index browser
def setup_roiURLs(get_roiURL, target_factors):
    global roiURLs
    roiURLs = []
    for dl in range(1, len(included_factors) + 1):
        if get_roiURL == 'Tru':  # 'True'
            # roiURLs.append('roiURL not available')
            roiURLs.append(roiURL(roi, top_region_coords[dl-1],
                                  dl, '_'.join(target_factors.astype(str))))
        else:
            roiURLs.append('roiURL not available')


# reorder data based on ranks in query region/roi
def reorder_by_query_region():
    global return_top_regions_per, return_pval, return_Z
    return_top_regions_per = top_regions_per[:, :, query_factor_order]
    return_pval = pval[:, query_factor_order]
    return_Z = Z[:, query_factor_order]


# build array of objects for roi and each similar region at each detail level
def generate_result_array(num_similar_regions):
    global return_arr_full
    return_arr_full = []
    for dl in range(len(included_factors)):
        # for the roi, return results as json array
        # each element is an object of scores for each metric
        return_arr = []
        for rank in range(num_similar_regions + 1):
            coord = top_region_coords[dl, rank]
            region_dict = {
                'rank': rank,
                'coordinates': coord[0] + ':' + coord[1] + '-' + coord[2],
                'region_ind': int(top_regions[dl, rank]),
                'score': float(top_regions_sim[dl, rank]),
                'percentiles': return_top_regions_per[dl, rank].tolist(),
                'prod_rank': query_sim_prod_ranks[dl, rank]
            }
            if rank == 0:
                region_dict['region_ind'] = 'NA'
                region_dict['score'] = 'NA'
                region_dict['factor_order'] = query_factor_order.tolist()
                region_dict['target_factors'] = target_factors.tolist()
                region_dict['top_region_factor_Z_score'] = return_Z[dl].tolist()
                region_dict['top_region_factor_pval'] = return_pval[dl].tolist()
                region_dict['roiURL'] = roiURLs[dl]
                region_dict['detail_level'] = dl + 1
                # region_dict['query_prod_rank'] = query_sim_prod_ranks[dl]
            return_arr.append(region_dict)

        return_arr_full.append(return_arr)


# find the factors we want to search over
def setup_selected_factors(selected_factors_str):
    global selected_factors
    selected_factors = np.int_(
        selected_factors_str.strip('[').strip(']').split(',')
    )


def determine_factor_order(top_regions, selected_factors):
    top_region_ranks = percents[top_regions]

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
            percents[:, selected_factors_best_order],
            axis=1
        ) > 0
    )[0]

    total_factors = np.arange(percents.shape[1])
    leftover_factors = total_factors[~np.isin(total_factors, selected_factors_best_order)]
    leftover_factors_best_order = leftover_factors[
        np.argsort(
            -np.sum(
                percents[regions_to_consider][:, leftover_factors],
                axis=0
            ),
            kind='stable'
        )
    ]

    factor_order = np.append(selected_factors_best_order, leftover_factors_best_order)
    return factor_order


# search for the best regions for selected factors
def search_over_selected_factors(num_regions):
    global top_regions, top_regions_ranks, top_regions_coords, top_regions_score, factor_order
    # focus on selected factors
    selected_factor_search_space_ranks = percents[:, selected_factors]
    # fusion
    selected_factor_search_space_score = 1 - np.prod(1 - selected_factor_search_space_ranks, axis=1)
    # get the best regions
    top_regions = np.argsort(-selected_factor_search_space_score)[:num_regions]

    # factor order
    factor_order = determine_factor_order(top_regions, selected_factors)

    # find info for each top region
    top_regions_ranks = percents[top_regions][:, factor_order]
    top_regions_coords = coords[top_regions]
    top_regions_score = selected_factor_search_space_score[top_regions]


# build array of objects for each region at each detail level
def generate_search_by_factor_result_array(num_regions):
    global search_by_factor_return_arr
    # for each region, return relevant information
    search_by_factor_return_arr = []
    for rank in range(num_regions):
        coord = top_regions_coords[rank]
        region_dict = {
            'rank': rank,
            'coordinates': coord[0] + ':' + coord[1] + '-' + coord[2],
            'region_ind': int(top_regions[rank]),
            'score': float(top_regions_score[rank]),
            'percentiles': top_regions_ranks[rank].tolist()
        }
        if rank == 0:
            region_dict['factor_order'] = factor_order.tolist()
            region_dict['selected_factors'] = selected_factors.tolist()
        search_by_factor_return_arr.append(region_dict)


@app.route('/', methods=['GET'])
def about():
    python_version = sys.version
    return "<html><body><div>region-narration ({})</div><div>Python {}</div></body></html>".format(RN_VERSION_STR, python_version)


@app.route('/narration', methods=['POST'])
def narration():
    request_info = request.get_json()
    genomic_location = request_info['location']
    scale = request_info['scale']
    included_factors_str = request_info['factors']
    get_roiURL = request_info['roiURL']
    query_factor_thresh = float(request_info['queryFactorThresh'])
    num_similar_regions = int(request_info['numSimilarRegions'])
    region_method = request_info['regionMethod']

    # get roi
    load_roi(genomic_location, build_tabix_path(scale, region_method))
    print(build_tabix_path(scale, region_method), flush=True)

    # choose search space
    choose_search_space(scale, region_method)

    # factors the user cares about
    setup_included_factors(included_factors_str, query_factor_thresh)

    # find most similar region inds
    find_most_similar_region_indices(num_similar_regions)

    # dual track roiURL
    setup_roiURLs(get_roiURL, target_factors)

    # reorder based on query region
    reorder_by_query_region()

    # generate result array
    generate_result_array(num_similar_regions)

    # result array as json
    json_return = json.dumps(return_arr_full)

    return json_return


# find the top regions for selected factors
@app.route('/search_by_factor', methods=['POST'])
def search_by_factor():
    request_info = request.get_json()
    scale = request_info['scale']
    selected_factors_str = request_info['factors']
    num_regions = int(request_info['numRegions'])
    region_method = request_info['regionMethod']

    # choose search space
    choose_search_space(scale, region_method)

    # get the selected factors
    setup_selected_factors(selected_factors_str)

    # get the top regions for the selected factors
    search_over_selected_factors(num_regions)

    # build return array
    generate_search_by_factor_result_array(num_regions)

    # return results
    json_return = json.dumps(search_by_factor_return_arr)
    return json_return


if __name__ == '__main__':
    app.run(host='0.0.0.0')
