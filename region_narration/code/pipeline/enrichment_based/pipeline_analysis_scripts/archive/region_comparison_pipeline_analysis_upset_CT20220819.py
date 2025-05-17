import numpy as np
import sys
import os
import time
from datetime import date
import h5py
import tabix
from scipy.stats import rankdata
import pickle
from region_comparison_pipeline_analysis_functions_multiple_CT20220818 import combine_scores_multiple, region_return_multiple, load_hdf5_files_multiple, load_hdf5_genes_multiple, create_arr_multiple, combine_scores_multiple_total
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import upsetplot
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
today = date.today()

def make_full_upset(top_region_inds, DHS_arr_top_100, TF_arr_top_100, 
                    epilogos_arr_top_100, string_arr_top_100, DHS_ooe, TF_ooe, 
                    epilogos_ooe, DHS_roi_rank, TF_roi_rank, epilogos_roi_mean_rank, 
                    epilogos_roi_max_rank, DHS_component_order, TF_component_order, 
                    epilogos_component_order, string_comp_order, DHSfact2descr, 
                    TFfact2descr, epilogosfact2descr, geneID2name, all_arr_ori, colors, 
                    regions_path, roi):
    
    color_tracker = np.array([])
    topfactor2region = {}
    #find top comps for each metric/recommended region
    for i in range(len(DHS_arr_top_100)):
        if DHS_arr_top_100[i] > 0:
            DHS_i_rank = rankdata(-DHS_ooe[top_region_inds[i+1]])
            top_comp = np.argsort((DHS_roi_rank+DHS_i_rank)/2)[0]
            f = DHSfact2descr[DHS_component_order[top_comp]]
            if f not in topfactor2region:
                topfactor2region[f] = np.array([i])
                color_tracker = np.append(color_tracker, colors[0])
            else:
                topfactor2region[f] = np.append(topfactor2region[f], [i])
    
    for i in range(len(TF_arr_top_100)):
        if TF_arr_top_100[i] > 0:
            TF_i_rank = rankdata(-TF_ooe[top_region_inds[i+1]])
            top_comp = np.argsort((TF_roi_rank+TF_i_rank)/2)[0]
            f = TF_component_order[top_comp] + ': ' + TFfact2descr[TF_component_order[top_comp]]
            if f not in topfactor2region:
                topfactor2region[f] = np.array([i])
                color_tracker = np.append(color_tracker, colors[1])
            else:
                topfactor2region[f] = np.append(topfactor2region[f], [i])
    
    for i in range(len(epilogos_arr_top_100)):
        if epilogos_arr_top_100[i] > 0:
            epilogos_i_mean_rank = rankdata(-epilogos_ooe[top_region_inds[i+1]][:18])
            epilogos_i_max_rank = rankdata(-epilogos_ooe[top_region_inds[i+1]][18:])
            top_mean_comp = np.argsort((epilogos_roi_mean_rank + epilogos_i_mean_rank)/2)[0]
            top_max_comp = np.argsort((epilogos_roi_max_rank + epilogos_i_max_rank)/2)[0]
            if top_mean_comp==top_max_comp:
                f = epilogosfact2descr[epilogos_component_order[top_mean_comp]]
            else:
                f = (epilogosfact2descr[epilogos_component_order[top_mean_comp]] + ', with ' + 
                     epilogosfact2descr[epilogos_component_order[top_max_comp]] + ' peaks')
            if f not in topfactor2region:
                topfactor2region[f] = np.array([i])
                color_tracker = np.append(color_tracker, colors[2])
            else:
                topfactor2region[f] = np.append(topfactor2region[f], [i])
    
    #STRING Scores
    for i in range(len(string_arr_top_100[0, :])):
        if np.sum(string_arr_top_100[:, i]) > 0:
            pos_string_inds = np.where(string_arr_top_100[:, i]>0)[0]
            for comp in pos_string_inds:
                ori_gene = geneID2name[all_arr_ori[comp, i]]
                f = 'Gene ' + string_comp_order[comp] + ' ' + ori_gene
                if f not in topfactor2region:
                    topfactor2region[f] = np.array([i])
                    color_tracker = np.append(color_tracker, colors[comp+3])
                else:
                    topfactor2region[f] = np.append(topfactor2region[f], [i])
                
    upset_data = upsetplot.from_contents(topfactor2region)
    upset_data = upset_data.groupby(list(upset_data.index.names)).count()['id']

    #sorting upset categories
    upset_factors = np.array(upset_data.index.names)
    upset_factor_counts = np.array([])
    for i in range(len(upset_factors)):
        upset_factor_counts = np.append(upset_factor_counts, np.sum(upset_data.unstack(level=i)[True]))
    upset_factors_sorted = upset_factors[np.argsort(upset_factor_counts)]
    sorted_upset_data = upset_data.reorder_levels(upset_factors_sorted)
    color_tracker_sorted = color_tracker[np.argsort(upset_factor_counts)]
    
    fig = plt.figure(figsize=(20,10))
    u = upsetplot.plot(sorted_upset_data, 
                       sort_categories_by=None, 
                       show_counts=True, 
                       min_subset_size=1, 
                       sort_by='cardinality', 
                       element_size=50)

    ###setting label colors to match bar plots
    for label, c in zip(u['matrix'].yaxis.get_ticklabels(), color_tracker_sorted):
        label.set_backgroundcolor(c)

    f = open(regions_path, 'r')
    regions = f.readlines()
    regions = regions[roi].strip().split('\t')
    genome_coords = regions[0] + ':' + regions[1] + '-' + regions[2]
    f.close()

    plt.ylabel('Number of Regions')
    title='Breakdown of Factors with Top 100 Regions for ROI %s (%s)' %(roi, genome_coords)
    plt.suptitle(title, fontsize=18)
    print(title)
    plt.savefig('/home/ctrader/region_similarity/data/figures/' + title.replace(' ', '_').replace(':', '').replace('(', '').replace(')', '').replace('\n', '_') + '_CT' + today.strftime("%Y%m%d") + '.pdf', 
                bbox_inches='tight', 
                dpi=400, 
                facecolor="#FFFFFF", 
                edgecolor="#FFFFFF", 
                transparent=False)
    return sorted_upset_data, topfactor2region, color_tracker_sorted

def make_limited_upset(limiting_per, sorted_upset_data, topfactor2region, 
                       color_tracker_sorted, regions_path, roi):
    overall_set = set()
    topfactor2region_limited = {}
    ind = 0
    while len(overall_set) < limiting_per:
        factor = np.flip(list(sorted_upset_data.index.names))[ind]
        factor_set = set(list(topfactor2region[factor]))
        overall_set.update(factor_set)
        ind += 1
        topfactor2region_limited[factor] = topfactor2region[factor]
    
    upset_data_limited = upsetplot.from_contents(topfactor2region_limited)
    upset_data_limited = upset_data_limited.groupby(list(upset_data_limited.index.names)).count()['id']
    
    if len(topfactor2region_limited) > 1: #more than 1 factor
        #sorting upset categories
        upset_factors_limited = np.array(upset_data_limited.index.names)
        upset_factor_counts_limited = np.array([])

        for i in range(len(upset_factors_limited)):
            upset_factor_counts_limited = np.append(upset_factor_counts_limited, np.sum(upset_data_limited.unstack(level=i)[True]))
        upset_factors_sorted_limited = upset_factors_limited[np.argsort(upset_factor_counts_limited)]
        sorted_upset_data_limited = upset_data_limited.reorder_levels(upset_factors_sorted_limited)

        fig = plt.figure(figsize=(20,10))

        u = upsetplot.plot(sorted_upset_data_limited, 
                           sort_categories_by=None, 
                           show_counts=True, 
                           min_subset_size=1, 
                           sort_by='cardinality', 
                           element_size=50)

        ###setting label colors to match bar plots
        for label, c in zip(u['matrix'].yaxis.get_ticklabels(), color_tracker_sorted[-len(upset_factors_sorted_limited):]):
            label.set_backgroundcolor(c)

        f = open(regions_path, 'r')
        regions = f.readlines()
        regions = regions[roi].strip().split('\t')
        genome_coords = regions[0] + ':' + regions[1] + '-' + regions[2]
        f.close()

        plt.ylabel('Number of Regions')
        title='Factors Describing ' + str(limiting_per) + '%' + ' of Top 100 Regions\nfor ROI %s (%s)' %(roi, genome_coords)
        plt.suptitle(title, fontsize=18)
        print(title)
        plt.savefig('/home/ctrader/region_similarity/data/figures/' + title.replace(' ', '_').replace(':', '').replace('(', '').replace(')', '').replace('\n', '_') + '_CT' + today.strftime("%Y%m%d") + '.pdf', 
                    bbox_inches='tight', 
                    dpi=400, 
                    facecolor="#FFFFFF", 
                    edgecolor="#FFFFFF", 
                    transparent=False)
    return


def main(gz_file_path, genomic_location, DHS_hdf5_path, TF_hdf5_path, epilogos_hdf5_path, coex_hdf5_path, cooccur_hdf5_path, database_hdf5_path, exper_hdf5_path,
         fusion_hdf5_path, homology_hdf5_path, neigh_hdf5_path, textmin_hdf5_path, regions_path, DHS_ooe_path, TF_ooe_path, epilogos_ooe_path,
         DHS_mapped_path, TF_mapped_path, epilogos_mapped_path, DHSfact2descr_path,  TFfact2descr_path, epifact2descr_path, geneID2name_path, out_bed_path, limiting_per):

    #load ooe files
    DHS_ooe_hdf5 = h5py.File(DHS_ooe_path, 'r')
    DHS_ooe = DHS_ooe_hdf5['ooe_arr'][:]
    DHS_ooe_hdf5.close()
    TF_ooe_hdf5 = h5py.File(TF_ooe_path, 'r')
    TF_ooe = TF_ooe_hdf5['ooe_arr'][:]
    TF_ooe_hdf5.close()
    epilogos_ooe_hdf5 = h5py.File(epilogos_ooe_path, 'r')
    epilogos_ooe = epilogos_ooe_hdf5['ooe_arr'][:]
    epilogos_ooe_hdf5.close()

    #descriptors for each factor
    pickle_in = open(TFfact2descr_path, 'rb')
    TFfact2descr = pickle.load(pickle_in)
    pickle_in = open(DHSfact2descr_path, 'rb')
    DHSfact2descr = pickle.load(pickle_in)
    pickle_in = open(epifact2descr_path, 'rb')
    epilogosfact2descr = pickle.load(pickle_in)

    #loading order of each component
    f = open(TF_mapped_path, 'r')
    TF_component_order = np.array(f.readline().strip().split('\t'))[3:]
    f.close()
    f = open(DHS_mapped_path, 'r')
    DHS_component_order = np.array(f.readline().strip().split('\t'))[3:]
    f.close()
    f = open(epilogos_mapped_path, 'r')
    epilogos_component_order = np.array(f.readline().strip().split('\t'))[3:]
    f.close()

    pickle_in = open(geneID2name_path, 'rb')
    geneID2name = pickle.load(pickle_in)

    #retrieve score arrays for regions of interest
    roi, top_region_inds, all_metrics_arr_top_100, all_arr_ori, all_arr_dest = combine_scores_multiple_total(gz_file_path, genomic_location, DHS_hdf5_path, TF_hdf5_path, epilogos_hdf5_path, 
                                                                                                             coex_hdf5_path, cooccur_hdf5_path, database_hdf5_path, exper_hdf5_path, fusion_hdf5_path, 
                                                                                                             homology_hdf5_path, neigh_hdf5_path, textmin_hdf5_path, regions_path, DHS_ooe_path, 
                                                                                                             TF_ooe_path, epilogos_ooe_path, DHS_mapped_path, TF_mapped_path, epilogos_mapped_path, 
                                                                                                             geneID2name_path, out_bed_path)
    #unpack
    score_arr_top_100 = all_metrics_arr_top_100[0]
    DHS_arr_top_100 = all_metrics_arr_top_100[1]
    TF_arr_top_100 = all_metrics_arr_top_100[2]
    epilogos_arr_top_100 = all_metrics_arr_top_100[3]
    string_arr_top_100 = all_metrics_arr_top_100[4:]

    #roi comps
    DHS_roi_rank = rankdata(-DHS_ooe[top_region_inds[:, 0]], axis=1)
    TF_roi_rank = rankdata(-TF_ooe[top_region_inds[:, 0]], axis=1)
    epilogos_roi_mean_rank = rankdata(-epilogos_ooe[top_region_inds[:, 0]][:, :18], axis=1)
    epilogos_roi_max_rank = rankdata(-epilogos_ooe[top_region_inds[:, 0]][:, 18:], axis=1)

    #sideways bar plots ordered by biggest component
    all_metrics_arr_top_100 = all_metrics_arr_top_100[1:]

    metric_sum_top_100 = np.sum(all_metrics_arr_top_100, axis=0)
    metrics_bar_scaled = all_metrics_arr_top_100 / metric_sum_top_100 * score_arr_top_100
    metrics_bar_scaled_inds = np.argsort(-metrics_bar_scaled, axis=0)
    metrics_bar_scaled_vals = -np.sort(-metrics_bar_scaled, axis=0)

    colors = np.array(['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', 
                       '#cab2d6', '#6a3d9a', '#ffff99'])
    string_comp_order = ['Co-expressed with', 'Co-occurs within Species with', 'Co-annotated with', 
                         'Interacts in vitro with', 'Fused in Other Species with', 'Has Common Ancestory with', 
                         'is Nearby', 'Co-occurs in Literature with']

    for r in range(len(roi)):
        (sorted_upset_data, 
         topfactor2region, 
         color_tracker_sorted) = make_full_upset(top_region_inds[r], DHS_arr_top_100[r], 
                                                 TF_arr_top_100[r], epilogos_arr_top_100[r], 
                                                 string_arr_top_100[:, r], DHS_ooe, TF_ooe, 
                                                 epilogos_ooe, DHS_roi_rank[r], TF_roi_rank[r], 
                                                 epilogos_roi_mean_rank[r], epilogos_roi_max_rank[r], 
                                                 DHS_component_order, TF_component_order, 
                                                 epilogos_component_order, string_comp_order, 
                                                 DHSfact2descr, TFfact2descr, epilogosfact2descr, 
                                                 geneID2name, all_arr_ori[:, r], colors, 
                                                 regions_path, roi[r])
        make_limited_upset(limiting_per, sorted_upset_data, topfactor2region, 
                           color_tracker_sorted, regions_path, roi[r])

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12],
         sys.argv[13], sys.argv[14], sys.argv[15], sys.argv[16], sys.argv[17], sys.argv[18], sys.argv[19], sys.argv[20], sys.argv[21], sys.argv[22], sys.argv[23], 
         sys.argv[24], sys.argv[25], int(sys.argv[26]))
