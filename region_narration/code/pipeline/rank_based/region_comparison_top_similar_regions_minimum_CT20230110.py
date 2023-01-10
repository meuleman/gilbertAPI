import numpy as np
import os
from datetime import date
import sys
import h5py
import tabix
import json
from region_comparison_sentences_CT20221129 import make_sentences


def combine_scores_multiple(DHS_scaled_vals, DHS_similar_regions,
                            TF_scaled_vals, TF_similar_regions,
                            epi_scaled_vals, epi_similar_regions):
    # max ind for combined score array
    max_ind = np.max([DHS_similar_regions, TF_similar_regions,
                      epi_similar_regions])

    # DHS scores
    DHS_arr = create_arr_multiple(DHS_similar_regions, DHS_scaled_vals,
                                  max_ind)
    # epilogos scores
    epi_arr = create_arr_multiple(epi_similar_regions, epi_scaled_vals,
                                  max_ind)
    # TF motif scores
    TF_arr = create_arr_multiple(TF_similar_regions, TF_scaled_vals, max_ind)

    # scoring all metrics
    chromatin_arr = 1 - ((1 - epi_arr) * (1 - DHS_arr) * (1 - TF_arr))

    # labeling most contributing metrics for each region
    contributing_chromatin_metrics_arr = np.argmax(np.array((DHS_arr, TF_arr,
                                                             epi_arr)),
                                                   axis=0)

    return (chromatin_arr, contributing_chromatin_metrics_arr,
            max_ind, DHS_arr, TF_arr, epi_arr)


def region_return_multiple(score_arr, roi):

    # top 100 region inds and scores
    rec_inds = np.argsort(-score_arr)[:, :100]
    top_region_inds = np.hstack((roi.reshape(-1, 1), rec_inds))
    # top_region_scores = np.append(1, score_arr[rec_inds])

    return top_region_inds


def roi_url(roi, top_regions):
    today = date.today()
    roi_url_name = ('region_comparison_recommended_regions_roi_' + str(roi) +
                    '_scale_25kb_CT' + today.strftime("%Y%m%d") + '.bed')
    bed_path_name = os.path.join('/home/ctrader/public_html/', roi_url_name)
    np.savetxt(bed_path_name, top_regions, delimiter='\t', fmt='%s')
    link_path1 = ('https://epilogos.altius.org:1234/?application=viewer&' +
                  'sampleSet=vC&mode=single&genome=hg38&model=18&complex' +
                  'ity=KL&group=all&chrLeft=chr1&chrRight=chr1&start=800' +
                  '0&stop=35000&roiURL=https%3A%2F%2Fencode%3Acollabor8%' +
                  '40resources.altius.org%2F~ctrader%2F' + roi_url_name)
    link_path2 = ('https://index.altius.org:1234/?application=viewer&mod' +
                  'e=pd&chrLeft=chr11&chrRight=chr11&start=47329800&stop' +
                  '=47356800&roiURL=https%3A%2F%2Fencode%3Acollabor8%40r' +
                  'esources.altius.org%2F~ctrader%2F' + roi_url_name)
    return link_path1, link_path2, bed_path_name


def load_hdf5_files_multiple(hdf5_path, roi):
    # take in a file path and load scaled scores and similar regions
    hdf5_file = h5py.File(hdf5_path, 'r')
    similar_regions = hdf5_file['top_1000_similar_regions'][roi]
    factor_ranks = hdf5_file['ranked_regions'][roi]
    sigmoid_scaled_vals = hdf5_file['top_1000_sigmoid_scaled_vals'][roi]
    hdf5_file.close()
    return sigmoid_scaled_vals, similar_regions, factor_ranks


def create_arr_multiple(similar_regions, scaled_vals, max_ind):
    # empty array
    arr = np.zeros((len(similar_regions), max_ind + 1))

    # find the values to change
    cols = similar_regions
    rows = np.array([np.arange(len(similar_regions))] *
                    len(similar_regions[0])).T

    # set values
    arr[rows, cols] = scaled_vals
    return arr


def get_factor_order(path):
    if path[-5:] == '.hdf5':
        hdf5_file = h5py.File(path, 'r')
        factor_order = np.char.decode(hdf5_file['comp_order'][:])
        hdf5_file.close()
    else:
        f = open(path, 'r')
        factor_order = np.array(f.readline().strip().split('\t'))[3:]
        f.close()
    return factor_order


def load_rank_path(hdf5_path, roi):
    hdf5_file = h5py.File(hdf5_path, 'r')
    rank_arr = hdf5_file['ranked_regions'][roi]
    max_ranks = hdf5_file['max_ranks'][:]
    hdf5_file.close()
    return rank_arr, max_ranks


def main(DHS_hdf5_path, TF_hdf5_path, epilogos_hdf5_path,
         DHS_mapped_path, TF_mapped_path, epilogos_mapped_path,
         STRING_ranks_hdf5_path, gz_file_path,
         regions_path, genomic_location, search_method):

    # load roi
    gz_tb = tabix.open(gz_file_path)
    roi = np.array(list(gz_tb.querys(genomic_location)))[:, 3].astype(int)
 
    # Load rank and ind arrays
    # DHS counts
    (DHS_scaled_vals, DHS_similar_regions,
     DHS_roi_ranks) = load_hdf5_files_multiple(DHS_hdf5_path, roi)
    # TF motif counts
    (TF_scaled_vals, TF_similar_regions,
     TF_roi_ranks) = load_hdf5_files_multiple(TF_hdf5_path, roi)
    # epilogos score sum and max
    (epi_scaled_vals, epi_similar_regions,
     epi_roi_ranks) = load_hdf5_files_multiple(epilogos_hdf5_path, roi)

    # Combine scores
    (chromatin_arr,
     contributing_chromatin_metrics_arr,
     max_ind, DHS_arr, TF_arr,
     epi_arr) = combine_scores_multiple(DHS_scaled_vals,
                                        DHS_similar_regions,
                                        TF_scaled_vals,
                                        TF_similar_regions,
                                        epi_scaled_vals,
                                        epi_similar_regions)

    # Return regions
    if search_method == 'chromatin':
        top_region_inds = region_return_multiple(chromatin_arr, roi)

    # values for top 100 regions (including query region)
    cols = top_region_inds
    rows = np.array([np.arange(len(top_region_inds))] *
                    len(top_region_inds[0])).T

    chromatin_arr_top_100 = np.array((DHS_arr, TF_arr,
                                      epi_arr))[:, rows, cols]
    chromatin_metric_order = ['DHS', 'TF', 'epilogos']

    # genomic coordinates for each region
    f = open(regions_path, 'r')
    regions = []
    for line in f.readlines():
        line = line.strip().split('\t')
        regions.append(line[0] + ':' + line[1] + '-' + line[2])
    f.close()

    top_region_coords = np.array(regions)[top_region_inds]

    sort_order = np.argsort(top_region_inds[0])
    # sentences for all regions
    sentences = make_sentences(DHS_hdf5_path, DHS_mapped_path, TF_hdf5_path,
                               TF_mapped_path, epilogos_hdf5_path,
                               epilogos_mapped_path, STRING_ranks_hdf5_path,
                               roi=top_region_inds[0][sort_order])
    sentences = sentences[np.argsort(sort_order)]

    # for the first roi, return results as json array
    # each element is an object of scores for each metric
    chromatin_return = []
    for rank in range(101):
        region_dict = {}
        for metric, name in zip(range(len(chromatin_metric_order)),
                                chromatin_metric_order):
            region_dict[name] = chromatin_arr_top_100[metric, 0, rank]
        region_dict['rank'] = rank
        region_dict['coordinates'] = top_region_coords[0, rank]
        region_dict['region_ind'] = int(top_region_inds[0, rank])
        region_dict['sentence'] = sentences[rank]
        chromatin_return.append(region_dict)

    json_return = json.dumps(chromatin_return)
    sys.stdout.write(json_return)

    # save output for dev mode
    # f = open('/Users/ctrader/Documents/AG/NarrativeGenomics/' +
    #          'region_narration/code/pipeline/region_similarity_react/' +
    #          'region-narration-app/src/' +
    #          'devmode_top_similar_regions_results_CT20221129.json', 'w')
    # f.write(json_return)
    # f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
         sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10],
         sys.argv[11], sys.argv[12], sys.argv[13], sys.argv[14], sys.argv[15],
         sys.argv[16], sys.argv[17], sys.argv[18], sys.argv[19])
