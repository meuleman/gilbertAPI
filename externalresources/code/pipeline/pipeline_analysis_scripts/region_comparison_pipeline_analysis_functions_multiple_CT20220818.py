import numpy as np
import sys
import os
import time
from datetime import date
import h5py
import tabix
from scipy.stats import rankdata
import pickle


def combine_scores_multiple(DHS_count_scaled_vals, DHS_count_similar_regions, TF_count_scaled_vals, TF_count_similar_regions, epilogos_score_scaled_vals,
                            epilogos_score_similar_regions, coex_scaled_vals, coex_similar_regions, cooccur_scaled_vals, cooccur_similar_regions,
                            database_scaled_vals, database_similar_regions, exper_scaled_vals, exper_similar_regions,fusion_scaled_vals, fusion_similar_regions,
                            homology_scaled_vals, homology_similar_regions, neigh_scaled_vals, neigh_similar_regions, textmin_scaled_vals, textmin_similar_regions):
    #max ind for combined score array
    max_ind = np.max([DHS_count_similar_regions, TF_count_similar_regions, epilogos_score_similar_regions, coex_similar_regions,
                     cooccur_similar_regions, database_similar_regions, exper_similar_regions, fusion_similar_regions, homology_similar_regions,
                     neigh_similar_regions, textmin_similar_regions])

    #DHS scores
    DHS_arr = create_arr_multiple(DHS_count_similar_regions, DHS_count_scaled_vals, max_ind)
    #epilogos scores
    epilogos_arr = create_arr_multiple(epilogos_score_similar_regions, epilogos_score_scaled_vals, max_ind)
    #TF motif scores
    TF_arr = create_arr_multiple(TF_count_similar_regions, TF_count_scaled_vals, max_ind)
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

    #scoring all metrics
    score_arr = 1 - ((1 - epilogos_arr) * (1 - DHS_arr) * (1 - TF_arr) * (1 - coex_arr) * (1 - cooccur_arr) * (1 - database_arr) * (1 - exper_arr) *
                     (1 - fusion_arr) * (1 - homology_arr) * (1 - neigh_arr) * (1 - textmin_arr))

    #labeling most contributing metrics for each region
    contributing_cosine_arr = np.argmax(np.array((DHS_arr, TF_arr, epilogos_arr, coex_arr, cooccur_arr, database_arr, exper_arr, fusion_arr,
                                                  homology_arr, neigh_arr, textmin_arr)), axis=0)

    return (score_arr, contributing_cosine_arr, max_ind, DHS_arr, TF_arr, epilogos_arr, coex_arr, 
            cooccur_arr, database_arr, exper_arr, fusion_arr, homology_arr, neigh_arr, textmin_arr)

def region_return_multiple(score_arr, contributing_cosine_arr, roi, regions_path, DHS_ooe_path, TF_ooe_path, epilogos_ooe_path, DHS_mapped_path, TF_mapped_path,
                           epilogos_mapped_path, coex_genes, cooccur_genes, database_genes, exper_genes, fusion_genes, homology_genes, neigh_genes,
                           textmin_genes, geneID2name_path):

    #top 100 region inds and scores
    rec_inds = np.argsort(-score_arr)[:, :100]
    top_region_inds = np.hstack((roi.reshape(-1, 1), rec_inds))
    # top_region_scores = np.append(1, score_arr[rec_inds])

    return top_region_inds

def load_hdf5_files_multiple(hdf5_path, roi):
    #take in a file path and load scaled scores and similar regions
    hdf5_file = h5py.File(hdf5_path, 'r')
    scaled_vals = hdf5_file['top_1000_scaled_vals'][roi]
    similar_regions = hdf5_file['top_1000_similar_regions'][roi]
    hdf5_file.close()
    return scaled_vals, similar_regions

def load_hdf5_genes_multiple(hdf5_path, roi, similar_regions, max_ind):
    #load in tracking genes
    hdf5_file = h5py.File(hdf5_path, 'r')
    ori_genes = np.char.decode(hdf5_file['top_1000_original_genes'][roi], 'utf-8')
    dest_genes = np.char.decode(hdf5_file['top_1000_destination_genes'][roi], 'utf-8')
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

def create_arr_multiple(similar_regions, scaled_vals, max_ind):
    #empty array
    arr = np.zeros((len(similar_regions), max_ind + 1))
    
    #find the values to change
    cols = similar_regions
    rows = np.array([np.arange(len(similar_regions))]*len(similar_regions[0])).T
    
    #set values
    arr[rows, cols] = scaled_vals
    return arr

def combine_scores_multiple_total(gz_file_path, genomic_location, DHS_hdf5_path, TF_hdf5_path, epilogos_hdf5_path, coex_hdf5_path, cooccur_hdf5_path, database_hdf5_path, exper_hdf5_path,
                                  fusion_hdf5_path, homology_hdf5_path, neigh_hdf5_path, textmin_hdf5_path, regions_path, DHS_ooe_path, TF_ooe_path, epilogos_ooe_path,
                                  DHS_mapped_path, TF_mapped_path, epilogos_mapped_path, geneID2name_path, out_bed_path):

    #find regions of interest    
    gz_tb = tabix.open(gz_file_path)
    roi = np.array(list(gz_tb.querys(genomic_location)))[:, 3].astype(int)
    print('roi = ' + str(roi), flush=True)

    ###Load rank and ind arrays
    t = time.time()

    #DHS counts
    DHS_count_scaled_vals, DHS_count_similar_regions = load_hdf5_files_multiple(DHS_hdf5_path, roi)
    #TF motif counts
    TF_count_scaled_vals, TF_count_similar_regions = load_hdf5_files_multiple(TF_hdf5_path, roi)
    #epilogos score sum and max
    epilogos_score_scaled_vals, epilogos_score_similar_regions = load_hdf5_files_multiple(epilogos_hdf5_path, roi)
    #coexpression STRING scores
    coex_scaled_vals, coex_similar_regions = load_hdf5_files_multiple(coex_hdf5_path, roi)
    #cooccurence STRING scores
    cooccur_scaled_vals, cooccur_similar_regions = load_hdf5_files_multiple(cooccur_hdf5_path, roi)
    #database STRING scores
    database_scaled_vals, database_similar_regions = load_hdf5_files_multiple(database_hdf5_path, roi)
    #experiments STRING scores
    exper_scaled_vals, exper_similar_regions = load_hdf5_files_multiple(exper_hdf5_path, roi)
    #fusion STRING scores
    fusion_scaled_vals, fusion_similar_regions = load_hdf5_files_multiple(fusion_hdf5_path, roi)
    #homology STRING scores
    homology_scaled_vals, homology_similar_regions = load_hdf5_files_multiple(homology_hdf5_path, roi)
    #neighborhood STRING scores
    neigh_scaled_vals, neigh_similar_regions = load_hdf5_files_multiple(neigh_hdf5_path, roi)
    #textmining STRING scores
    textmin_scaled_vals, textmin_similar_regions = load_hdf5_files_multiple(textmin_hdf5_path, roi)

    ###Combine scores 
    (score_arr, contributing_cosine_arr,
     max_ind, DHS_arr, TF_arr, epilogos_arr,
     coex_arr, cooccur_arr, database_arr,
     exper_arr, fusion_arr, homology_arr,
     neigh_arr, textmin_arr) = combine_scores_multiple(DHS_count_scaled_vals, DHS_count_similar_regions,
                                                       TF_count_scaled_vals, TF_count_similar_regions,
                                                       epilogos_score_scaled_vals, epilogos_score_similar_regions,
                                                       coex_scaled_vals, coex_similar_regions,
                                                       cooccur_scaled_vals, cooccur_similar_regions,
                                                       database_scaled_vals, database_similar_regions,
                                                       exper_scaled_vals, exper_similar_regions,
                                                       fusion_scaled_vals, fusion_similar_regions,
                                                       homology_scaled_vals, homology_similar_regions,
                                                       neigh_scaled_vals, neigh_similar_regions,
                                                       textmin_scaled_vals, textmin_similar_regions)

    #load in gene names for STRING scores
    coex_genes = load_hdf5_genes_multiple(coex_hdf5_path, roi, coex_similar_regions, max_ind)
    cooccur_genes = load_hdf5_genes_multiple(cooccur_hdf5_path, roi, cooccur_similar_regions, max_ind)
    database_genes = load_hdf5_genes_multiple(database_hdf5_path, roi, database_similar_regions, max_ind)
    exper_genes = load_hdf5_genes_multiple(exper_hdf5_path, roi, exper_similar_regions, max_ind)
    fusion_genes = load_hdf5_genes_multiple(fusion_hdf5_path, roi, fusion_similar_regions, max_ind)
    homology_genes = load_hdf5_genes_multiple(homology_hdf5_path, roi, homology_similar_regions, max_ind)
    neigh_genes = load_hdf5_genes_multiple(neigh_hdf5_path, roi, neigh_similar_regions, max_ind)
    textmin_genes = load_hdf5_genes_multiple(textmin_hdf5_path, roi, textmin_similar_regions, max_ind)

    ###Return regions
    top_region_inds = region_return_multiple(score_arr, contributing_cosine_arr, roi, regions_path, DHS_ooe_path, TF_ooe_path, epilogos_ooe_path, DHS_mapped_path, TF_mapped_path,
                                             epilogos_mapped_path, coex_genes, cooccur_genes, database_genes, exper_genes, fusion_genes, homology_genes, neigh_genes,
                                             textmin_genes, geneID2name_path)

    #taking the top 100 regions from each score array

    cols = top_region_inds[:, 1:]
    rows = np.array([np.arange(len(top_region_inds[:, 1:]))]*len(top_region_inds[:, 1:][0])).T

    all_metrics_arr_top_100 = np.array((score_arr, DHS_arr, TF_arr, epilogos_arr, coex_arr, 
                                        cooccur_arr, database_arr, exper_arr, fusion_arr, 
                                        homology_arr, neigh_arr, textmin_arr))[:, rows, cols]

    #all string original and destination genes
    all_arr_ori = np.array([coex_genes[0], cooccur_genes[0], database_genes[0], exper_genes[0], 
                            fusion_genes[0], homology_genes[0], neigh_genes[0], textmin_genes[0]])[:, rows, cols]
    all_arr_dest = np.array([coex_genes[1], cooccur_genes[1], database_genes[1], exper_genes[1], 
                             fusion_genes[1], homology_genes[1], neigh_genes[1], textmin_genes[1]])[:, rows, cols]

    print('Time to load data and find similar regions = ' + str(time.time()-t), flush=True)

    return roi, top_region_inds, all_metrics_arr_top_100, all_arr_ori, all_arr_dest
