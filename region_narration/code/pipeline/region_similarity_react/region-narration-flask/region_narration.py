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
prod_name = '/data/narrativeGenomics/NarrativeGenomics/region_narration/code/pipeline/region_similarity_react/region-narration-flask/region_narration.py'

if current_name == prod_name:
    server_root_dir = '/data/narrativeGenomics/NarrativeGenomics/region_narration'

    RN_VERSION = '0.0.1'
    RN_VERSION_STR = "{}/{}".format(Repo(server_root_dir, search_parent_directories=True).head.ref.name, RN_VERSION)

    # load search spaces for each scale
    data_dir = os.path.join(server_root_dir, 'data')

else:
    RN_VERSION_STR = ''
    # load search spaces for each scale
    data_dir = '../../../../data/'

# Given a scale, build a path name for the search space
def build_search_space_path(scale):
    path = os.path.join(data_dir, 'regions_{}kb_CT20230427.DHS.factornormNMF.sum_enrich_corank.epilogos.statesum_rank.percentiles.limited_search_space_top_1000_no_overlap.npz'.format(scale))
    return path


"""
# Given a scale, build a path name for the query region
def build_query_path(scale):
    path = os.path.join(data_dir, 'regions_{}kb_CT20230412.DHS.factornormNMF.sum_enrich_corank.epilogos.statesum_rank.percentiles.hdf5'.format(scale))
    return path
"""


# Given a scale, build a path name for the tabix query
def build_tabix_path(scale):
    # path = os.path.join(data_dir, 'hg38.{}kb.regions.CT20230412.sorted.rownum4.bed.gz'.format(scale))
    path = 'https://resources.altius.org/~ctrader/public/genomic_narration/indexed/hg38.{}kb.regions.CT20230427.sorted.rownum4.factorranks5t38.bed.gz'.format(scale)
    return path


# Load search space data from a given scale
def load_ranks(search_space_path):
    search_space = np.load(search_space_path)
    percents = search_space['percentiles']
    coords = search_space['coordinates']
    return percents, coords


# load search spaces for all scales into dict
seach_spaces = {}
for scale in [1, 2, 5, 10, 15, 25, 50, 75, 100]:
    seach_spaces[scale] = load_ranks(build_search_space_path(scale))


# Functions for narration

def tabix_query(filename, query):
    # Use subprocess with tabix and return hits
    process = Popen(['tabix', '-f', filename, query], stdout=PIPE)

    hits = []
    for line in process.stdout:
        hits.append(np.char.decode(line.strip().split()))
    hits = np.array(hits)
    return hits


"""
# determine roi from query coordinates
def load_roi(genomic_location, tabix_file_path):
    global roi

    # determine chromosome, start, and end from query
    query_start_end = genomic_location.split(':')[1].split('-')
    query_start = int(query_start_end[0])
    query_end = int(query_start_end[1])

    # retrieve possible rois
    tb = tabix.open(tabix_file_path)
    possible_coords = np.array(list(tb.querys(genomic_location)))

    # find closest roi
    start_dist = np.abs(possible_coords[:, 1].astype(int) - query_start)
    end_dist = np.abs(possible_coords[:, 2].astype(int) - query_end)
    query_region = possible_coords[np.argmin(start_dist + end_dist)]
    roi = int(query_region[3])

    return
"""


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
    query_region = possible_coords[np.argmin(start_dist + end_dist)]

    roi = int(query_region[3])
    roi_percents = query_region[4:].astype('float32')
    roi_coords = query_region[:3]

    return


"""
# load query ranks and coordinates
def load_ranks(query_path):
    global roi_percents, roi_coords

    # load roi ranks
    hdf5_f = h5py.File(query_path, 'r')
    roi_percents = hdf5_f['percentiles'][roi]
    roi_coords = np.char.decode(hdf5_f['coordinates'][roi])
    hdf5_f.close()
"""


# choose search space from scale options
def choose_search_space(scale):
    global percents, coords

    percents, coords = seach_spaces[scale]


# find the factors we want to search over
def setup_included_factors(included_factors_str, query_factor_thresh):
    global included_factors

    included_factors_from_str = np.int_(
        included_factors_str.strip('[').strip(']').split(',')
    )
    query_factors_above_thresh = np.where(roi_percents > query_factor_thresh)[0]

    # if no factors above threshold, search best one
    if len(query_factors_above_thresh) == 0:
        sort_order = np.argsort(-roi_percents)
        included_factors = sort_order[
            np.isin(sort_order, included_factors_from_str)
        ][:1]
    else:
        included_factors = included_factors_from_str[
            np.isin(included_factors_from_str, query_factors_above_thresh)
        ]


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
                region_dict['query_factor_order'] = query_factor_order.tolist()
                region_dict['target_factors'] = target_factors.tolist()
                region_dict['top_region_factor_Z_score'] = return_Z[dl].tolist()
                region_dict['top_region_factor_pval'] = return_pval[dl].tolist()
                region_dict['roiURL'] = roiURLs[dl]
                region_dict['detail_level'] = dl + 1
                # region_dict['query_prod_rank'] = query_sim_prod_ranks[dl]
            return_arr.append(region_dict)

        return_arr_full.append(return_arr)


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

    # get roi
    load_roi(genomic_location, build_tabix_path(scale))

    # # get ranks and coordinates for the query region
    # load_ranks(build_query_path(scale))

    # choose search space
    choose_search_space(scale)

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


if __name__ == '__main__':
    app.run(host='0.0.0.0')
