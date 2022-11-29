import numpy as np
import os
from datetime import date
import sys
import h5py
import tabix
import json

def combine_scores_multiple(DHS_scaled_vals, DHS_similar_regions, TF_scaled_vals, TF_similar_regions, epi_scaled_vals,
                            epi_similar_regions, coex_scaled_vals, coex_similar_regions, cooccur_scaled_vals, cooccur_similar_regions,
                            database_scaled_vals, database_similar_regions, exper_scaled_vals, exper_similar_regions,fusion_scaled_vals, fusion_similar_regions,
                            homology_scaled_vals, homology_similar_regions, neigh_scaled_vals, neigh_similar_regions, textmin_scaled_vals, textmin_similar_regions):
    #max ind for combined score array
    max_ind = np.max([DHS_similar_regions, TF_similar_regions, epi_similar_regions, coex_similar_regions,
                     cooccur_similar_regions, database_similar_regions, exper_similar_regions, fusion_similar_regions, homology_similar_regions,
                     neigh_similar_regions, textmin_similar_regions])

    #DHS scores
    DHS_arr = create_arr_multiple(DHS_similar_regions, DHS_scaled_vals, max_ind)
    #epilogos scores
    epi_arr = create_arr_multiple(epi_similar_regions, epi_scaled_vals, max_ind)
    #TF motif scores
    TF_arr = create_arr_multiple(TF_similar_regions, TF_scaled_vals, max_ind)
    #coexpression scores
    coex_arr = create_arr_multiple(coex_similar_regions, coex_scaled_vals, max_ind)
    #cooccurrence scores
    cooccur_arr = create_arr_multiple(cooccur_similar_regions, cooccur_scaled_vals, max_ind)
    #database scores
    database_arr = create_arr_multiple(database_similar_regions, database_scaled_vals, max_ind)
    #experiment scores
    exper_arr = create_arr_multiple(exper_similar_regions, exper_scaled_vals, max_ind)
    #fusion scores
    fusion_arr = create_arr_multiple(fusion_similar_regions, fusion_scaled_vals, max_ind)
    #homology scores
    homology_arr = create_arr_multiple(homology_similar_regions, homology_scaled_vals, max_ind)
    #neighborhood scores
    neigh_arr = create_arr_multiple(neigh_similar_regions, neigh_scaled_vals, max_ind)
    #text mining scores
    textmin_arr = create_arr_multiple(textmin_similar_regions, textmin_scaled_vals, max_ind)

    #fuse STRING scores first
    string_arr = 1 - ((1 - coex_arr) * (1 - cooccur_arr) * (1 - database_arr) * (1 - exper_arr) *
                      (1 - fusion_arr) * (1 - homology_arr) * (1 - neigh_arr) * (1 - textmin_arr))

    #scoring all metrics
    chromatin_arr = 1 - ((1 - epi_arr) * (1 - DHS_arr) * (1 - TF_arr))

    #labeling most contributing metrics for each region
    contributing_chromatin_metrics_arr = np.argmax(np.array((DHS_arr, TF_arr, epi_arr)), axis=0)
    contributing_string_metrics_arr = np.argmax(np.array((coex_arr, cooccur_arr, database_arr,
                                                          exper_arr, fusion_arr, homology_arr,
                                                          neigh_arr, textmin_arr)), axis=0)

    return (chromatin_arr, string_arr, contributing_chromatin_metrics_arr, contributing_string_metrics_arr, max_ind, DHS_arr,
            TF_arr, epi_arr, coex_arr, cooccur_arr, database_arr, exper_arr, fusion_arr, homology_arr,
            neigh_arr, textmin_arr)

def region_return_multiple(score_arr, roi):

    #top 100 region inds and scores
    rec_inds = np.argsort(-score_arr)[:, :100]
    top_region_inds = np.hstack((roi.reshape(-1, 1), rec_inds))
    # top_region_scores = np.append(1, score_arr[rec_inds])

    return top_region_inds

def roi_url(roi, top_regions):
    today = date.today()
    roi_url_name = 'region_comparison_recommended_regions_roi_'+str(roi)+'_scale_25kb_CT' + today.strftime("%Y%m%d") + '.bed'
    bed_path_name = os.path.join('/home/ctrader/public_html/', roi_url_name)
    np.savetxt(bed_path_name, top_regions, delimiter='\t', fmt='%s')
    link_path1 = 'https://epilogos.altius.org:1234/?application=viewer&sampleSet=vC&mode=single&genome=hg38&model=18&complexity=KL&group=all&chrLeft=chr1&chrRight=chr1&start=8000&stop=35000&roiURL=https%3A%2F%2Fencode%3Acollabor8%40resources.altius.org%2F~ctrader%2F' + roi_url_name
    link_path2 = 'https://index.altius.org:1234/?application=viewer&mode=pd&chrLeft=chr11&chrRight=chr11&start=47329800&stop=47356800&roiURL=https%3A%2F%2Fencode%3Acollabor8%40resources.altius.org%2F~ctrader%2F' + roi_url_name
    return link_path1, link_path2, bed_path_name

def load_hdf5_files_multiple(hdf5_path, roi):
    #take in a file path and load scaled scores and similar regions
    hdf5_file = h5py.File(hdf5_path, 'r')
    similar_regions = hdf5_file['top_1000_similar_regions'][roi]
    factor_ranks = hdf5_file['ranked_regions'][roi]
    sigmoid_scaled_vals = hdf5_file['top_1000_sigmoid_scaled_vals'][roi]
    hdf5_file.close()
    return sigmoid_scaled_vals, similar_regions, factor_ranks

def load_STRING_hdf5_files_multiple(hdf5_path, roi):
    #load in tracking genes
    hdf5_file = h5py.File(hdf5_path, 'r')
    ori_genes = np.char.decode(hdf5_file['top_1000_original_genes'][roi], 'utf-8')
    dest_genes = np.char.decode(hdf5_file['top_1000_destination_genes'][roi], 'utf-8')
    similar_regions = hdf5_file['top_1000_similar_regions'][roi]
    scaled_vals = hdf5_file['top_1000_scaled_vals'][roi]
    hdf5_file.close()
    return scaled_vals, similar_regions, ori_genes, dest_genes

def create_arr_multiple(similar_regions, scaled_vals, max_ind):
    #empty array
    arr = np.zeros((len(similar_regions), max_ind + 1))

    #find the values to change
    cols = similar_regions
    rows = np.array([np.arange(len(similar_regions))]*len(similar_regions[0])).T

    #set values
    arr[rows, cols] = scaled_vals
    return arr

def load_hdf5_genes_multiple(ori_genes, dest_genes, similar_regions, max_ind):
    #establishing gene arrays
    ori_arr = np.where(np.zeros((len(similar_regions), max_ind + 1))==0, None, None)
    dest_arr = np.where(np.zeros((len(similar_regions), max_ind + 1))==0, None, None)

    #filling gene arrays
    #find the values to change
    cols = similar_regions
    rows = np.array([np.arange(len(similar_regions))]*len(similar_regions[0])).T

    #set values
    ori_arr[rows, cols] = ori_genes
    dest_arr[rows, cols] = dest_genes
    return np.array([ori_arr, dest_arr])

def get_factor_order(path):
    if path[-5:]=='.hdf5':
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

def main(DHS_hdf5_path, TF_hdf5_path, epilogos_hdf5_path, coex_hdf5_path, cooccur_hdf5_path, database_hdf5_path, exper_hdf5_path, fusion_hdf5_path, homology_hdf5_path, neigh_hdf5_path,
         textmin_hdf5_path, gz_file_path, regions_path, genomic_location, search_method):

    #load roi
    gz_tb = tabix.open(gz_file_path)
    roi = np.array(list(gz_tb.querys(genomic_location)))[:, 3].astype(int)

    ###Load rank and ind arrays
    #DHS counts
    DHS_scaled_vals, DHS_similar_regions, DHS_roi_ranks = load_hdf5_files_multiple(DHS_hdf5_path, roi)
    #TF motif counts
    TF_scaled_vals, TF_similar_regions, TF_roi_ranks = load_hdf5_files_multiple(TF_hdf5_path, roi)
    #epilogos score sum and max
    epi_scaled_vals, epi_similar_regions, epi_roi_ranks = load_hdf5_files_multiple(epilogos_hdf5_path, roi)
    #coexpression STRING scores
    coex_scaled_vals, coex_similar_regions, coex_ori, coex_dest = load_STRING_hdf5_files_multiple(coex_hdf5_path, roi)
    #cooccurence STRING scores
    cooccur_scaled_vals, cooccur_similar_regions, cooccur_ori, cooccur_dest = load_STRING_hdf5_files_multiple(cooccur_hdf5_path, roi)
    #database STRING scores
    database_scaled_vals, database_similar_regions, database_ori, database_dest = load_STRING_hdf5_files_multiple(database_hdf5_path, roi)
    #experiments STRING scores
    exper_scaled_vals, exper_similar_regions, exper_ori, exper_dest = load_STRING_hdf5_files_multiple(exper_hdf5_path, roi)
    #fusion STRING scores
    fusion_scaled_vals, fusion_similar_regions, fusion_ori, fusion_dest = load_STRING_hdf5_files_multiple(fusion_hdf5_path, roi)
    #homology STRING scores
    homology_scaled_vals, homology_similar_regions, homology_ori, homology_dest = load_STRING_hdf5_files_multiple(homology_hdf5_path, roi)
    #neighborhood STRING scores
    neigh_scaled_vals, neigh_similar_regions, neigh_ori, neigh_dest = load_STRING_hdf5_files_multiple(neigh_hdf5_path, roi)
    #textmining STRING scores
    textmin_scaled_vals, textmin_similar_regions, textmin_ori, textmin_dest = load_STRING_hdf5_files_multiple(textmin_hdf5_path, roi)

    ###Combine scores
    (chromatin_arr, string_arr,
     contributing_chromatin_metrics_arr,
     contributing_string_metrics_arr,
     max_ind, DHS_arr, TF_arr,
     epi_arr, coex_arr, cooccur_arr,
     database_arr, exper_arr, fusion_arr,
     homology_arr, neigh_arr, textmin_arr) = combine_scores_multiple(DHS_scaled_vals, DHS_similar_regions,
                                                                     TF_scaled_vals, TF_similar_regions,
                                                                     epi_scaled_vals, epi_similar_regions,
                                                                     coex_scaled_vals, coex_similar_regions,
                                                                     cooccur_scaled_vals, cooccur_similar_regions,
                                                                     database_scaled_vals, database_similar_regions,
                                                                     exper_scaled_vals, exper_similar_regions,
                                                                     fusion_scaled_vals, fusion_similar_regions,
                                                                     homology_scaled_vals, homology_similar_regions,
                                                                     neigh_scaled_vals, neigh_similar_regions,
                                                                     textmin_scaled_vals, textmin_similar_regions)

    ##load in gene names for STRING scores
    #coex_genes = load_hdf5_genes_multiple(coex_ori, coex_dest, coex_similar_regions, max_ind)
    #cooccur_genes = load_hdf5_genes_multiple(cooccur_ori, cooccur_dest, cooccur_similar_regions, max_ind)
    #database_genes = load_hdf5_genes_multiple(database_ori, database_dest, database_similar_regions, max_ind)
    #exper_genes = load_hdf5_genes_multiple(exper_ori, exper_dest, exper_similar_regions, max_ind)
    #fusion_genes = load_hdf5_genes_multiple(fusion_ori, fusion_dest, fusion_similar_regions, max_ind)
    #homology_genes = load_hdf5_genes_multiple(homology_ori, homology_dest, homology_similar_regions, max_ind)
    #neigh_genes = load_hdf5_genes_multiple(neigh_ori, neigh_dest, neigh_similar_regions, max_ind)
    #textmin_genes = load_hdf5_genes_multiple(textmin_ori, textmin_dest, textmin_similar_regions, max_ind)

    ###Return regions
    if search_method=='chromatin':
        top_region_inds = region_return_multiple(chromatin_arr, roi)
    elif search_method=='string':
        top_region_inds = region_return_multiple(string_arr, roi)

    #values for top 100 regions
    cols = top_region_inds[:, 1:]
    rows = np.array([np.arange(len(top_region_inds[:, 1:]))]*len(top_region_inds[:, 1:][0])).T

    chromatin_arr_top_100 = np.array((DHS_arr, TF_arr, epi_arr))[:, rows, cols]
    chromatin_metric_order = ['DHS', 'TF motifs', 'epilogos']

    string_arr_top_100 = np.array((coex_arr, cooccur_arr, database_arr, exper_arr, fusion_arr, homology_arr,
                                   neigh_arr, textmin_arr))[:, rows, cols]
    string_metric_order = ['Co-expression', 'Co-occurence', 'Co-annotation',
                           'Interactions in vitro', 'Fusion in Other Species', 'Common Ancestry',
                           'being Nearby', 'Co-citation']

    #genomic coordinates for each region
    f = open(regions_path, 'r')
    regions = []
    for line in f.readlines():
        line = line.strip().split('\t')
        regions.append(line[0] + ':' + line[1] + '-' + line[2])
    f.close()
    top_region_coords = np.array(regions)[top_region_inds[:, 1:]]

    #for the first roi, return results as json array
    #each element is an object of scores for each metric
    chromatin_return = []
    for region in range(100):
        region_dict = {}
        for metric, name in zip(range(len(chromatin_arr_top_100[:, 0])), chromatin_metric_order):
            region_dict[name] = chromatin_arr_top_100[metric, 0, region]
        region_dict['region'] = region
        region_dict['coordinates'] = top_region_coords[0, region]
        chromatin_return.append(region_dict)

    json_return = json.dumps(chromatin_return)
    sys.stdout.write(json_return)
    #f = open('/Users/ctrader/Documents/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/devmode_top_similar_regions_results_CT20221123.json', 'w')
    #f.write(json_return)
    #f.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13], sys.argv[14], sys.argv[15])
