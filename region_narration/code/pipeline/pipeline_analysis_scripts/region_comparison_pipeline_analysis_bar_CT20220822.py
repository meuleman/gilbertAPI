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
from sklearn.cluster import MeanShift
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
today = date.today()

def make_barplot(top_region_inds, metrics_bar_scaled_vals, metrics_bar_scaled_inds, DHS_roi_rank, TF_roi_rank, 
                 epilogos_roi_mean_rank, epilogos_roi_max_rank, DHS_ooe, TF_ooe, epilogos_ooe, 
                 DHS_component_order, TF_component_order, epilogos_component_order, all_arr_ori, all_arr_dest, 
                 geneID2name, score_arr_top_100, DHSfact2descr, TFfact2descr, epilogosfact2descr, roi, colors, figure_out_dir):
    #plotting
    fig, axes = plt.subplots(ncols=2, figsize=(16, 20))
    fig.patch.set_facecolor('white')
    axes = axes.ravel()
    
    for i in range(len(metrics_bar_scaled_vals)):
        axes[0].barh(np.flip(np.arange(100)), metrics_bar_scaled_vals[i], left=np.sum(metrics_bar_scaled_vals[:i], axis=0), 
                     color=colors[metrics_bar_scaled_inds[i]], height=1, align='center')
        
    #label top comps
    top_factors = []
    for i in range(len(metrics_bar_scaled_vals[0])):
        start_val = 0
        for metric, val in zip(metrics_bar_scaled_inds[:, i], metrics_bar_scaled_vals[:, i]):
            if ((metric==0) & (val>0)):
                DHS_i_rank = rankdata(-DHS_ooe[top_region_inds[i+1]])
                top_comp = np.argsort((DHS_roi_rank+DHS_i_rank)/2)[0]
                if val / len(DHS_component_order[top_comp]) > 0.012:
                    axes[0].text(start_val + val/10, 
                             len(metrics_bar_scaled_vals[0])-1.35-i, 
                             DHS_component_order[top_comp], 
                             fontsize=9)
                    top_factors.append(DHS_component_order[top_comp])
            elif ((metric==1) & (val>0)):
                TF_i_rank = rankdata(-TF_ooe[top_region_inds[i+1]])
                top_comp = np.argsort((TF_roi_rank+TF_i_rank)/2)[0]
                if val / len(TF_component_order[top_comp]) > 0.012:
                    axes[0].text(start_val + val/10, 
                             len(metrics_bar_scaled_vals[0])-1.35-i, 
                             TF_component_order[top_comp], 
                             fontsize=9)
                    top_factors.append(TF_component_order[top_comp])
            elif ((metric==2) & (val>0)):
                epilogos_i_mean_rank = rankdata(-epilogos_ooe[top_region_inds[i+1]][:18])
                epilogos_i_max_rank = rankdata(-epilogos_ooe[top_region_inds[i+1]][18:])
                top_mean_comp = np.argsort((epilogos_roi_mean_rank + epilogos_i_mean_rank)/2)[0]
                top_max_comp = np.argsort((epilogos_roi_max_rank + epilogos_i_max_rank)/2)[0]
                if val / len(epilogos_component_order[top_mean_comp] + '/' + epilogos_component_order[top_max_comp]) > 0.012:
                    axes[0].text(start_val + val/10, 
                             len(metrics_bar_scaled_vals[0])-1.35-i, 
                             epilogos_component_order[top_mean_comp] + '/' + epilogos_component_order[top_max_comp], 
                             fontsize=9)
                    top_factors.append(epilogos_component_order[top_mean_comp] + '/' + epilogos_component_order[top_max_comp])
            elif ((metric>2) & (val>0)):
                metric_ori_gene = geneID2name[all_arr_ori[metric-3][i]]
                metric_dest_gene = geneID2name[all_arr_dest[metric-3][i]]
                if val / len(metric_ori_gene + ' to ' + metric_dest_gene) > 0.0135:
                    axes[0].text(start_val + val/10, 
                             len(metrics_bar_scaled_vals[0])-1.35-i, 
                             metric_ori_gene + ' to ' + metric_dest_gene, 
                             fontsize=9)
                    top_factors.append(metric_ori_gene + ' to X')
            start_val+=val
        axes[0].text(1.01, 
             len(metrics_bar_scaled_vals[0])-1.35-i, 
             np.round(score_arr_top_100[i], 4), 
             fontsize=9)
    axes[0].text(1.01, len(metrics_bar_scaled_vals[0]), 'Overall Score', fontsize=9)
    
    #look up table
    top_factors = np.unique(top_factors)
    factor_bar_arr = []
    for factor in top_factors:
        factor_l = factor.split('/')
        if len(factor_l) == 2: #epilogos
            if factor_l[0]==factor_l[1]:
                factor_bar_arr.append(epilogosfact2descr[factor_l[0]])
            else:
                factor_bar_arr.append(epilogosfact2descr[factor_l[0]] + ', with ' + epilogosfact2descr[factor_l[1]] + ' peaks')
        elif factor in DHSfact2descr: #DHS
            factor_bar_arr.append(DHSfact2descr[factor])
        elif factor in TFfact2descr: #TF motifs
            factor_bar_arr.append(factor + ': ' + TFfact2descr[factor])
        else:
            gene0, gene1 = factor.split(' to ')[0], factor.split(' to ')[1]
            factor_bar_arr.append('Gene ' + gene0 + ' Mapping to Gene X')
    factor_bar_arr = np.array(factor_bar_arr)
    lookup_table = np.hstack((top_factors.reshape((-1, 1)), factor_bar_arr.reshape((-1, 1))))
    t = axes[1].table(lookup_table, 
                      cellLoc='left',
                      colLabels=['Key', 'Descriptor'],
#                       bbox=[0, 0.3, 0.95, 0.7], 
                      bbox=[0.2, 1-(len(factor_bar_arr)/20), 0.75, (len(factor_bar_arr)/20)],
                      colWidths=[0.15, 0.35])
    t.auto_set_font_size(False)
    t.set_fontsize(9)

    axes[1].set_xticks([])
    axes[1].set_yticks([])

    axes[1].spines['right'].set_visible(False)
    axes[1].spines['left'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['bottom'].set_visible(False)

    axes[0].set_yticks(np.arange(9, 101, 10))
    axes[0].set_yticklabels(np.flip(np.arange(1, 92, 10)))
    axes[0].set_ylim([-10, 99.5])
    axes[0].set_xticks([1])
    h = []

    for metric, c in zip(['DNase I Accessibility', 'TF Motifs', 'Chromatin States', 'Co-expression', 
                          'Co-evolved', 'Co-annotated', 'Interacts', 'Fused', 'Homologous', 
                          'Neighborhood', 'Co-cited'], colors):
        p = mpatches.Patch(color=c, label=metric)
        h.append(p)
    axes[0].legend(handles=h, loc='lower center', ncol=3)
    axes[0].set_xlabel('Overall Score')
    axes[0].set_ylabel('Recommended Region')
    title = 'Scoring Breakdown for ROI %s' %(roi)
    axes[0].set_title(title)
    plt.savefig(figure_out_dir + title.replace(' ', '_') + '_CT' + today.strftime("%Y%m%d") + '.pdf', 
                bbox_inches='tight', 
                dpi=400, 
                facecolor="#FFFFFF", 
                edgecolor="#FFFFFF", 
                transparent=False)
    return

def make_cluster_barplot(all_metrics_arr_top_100, top_region_inds, score_arr_top_100, all_arr_ori, all_arr_dest, colors, 
                         DHS_roi_rank, TF_roi_rank, epilogos_roi_mean_rank, epilogos_roi_max_rank, DHS_ooe, 
                         TF_ooe, epilogos_ooe, DHS_component_order, TF_component_order, epilogos_component_order, 
                         geneID2name, DHSfact2descr, TFfact2descr, epilogosfact2descr, roi, figure_out_dir):
    #sideways bar plots (by clustering order)
    #cluster rows and get new order for regions
    cluster = MeanShift()
    cluster.fit(all_metrics_arr_top_100.T)
    cluster_labels = cluster.labels_
    new_order_top_100 = np.array([], dtype=int)
    for c in np.unique(cluster_labels):
        new_order_top_100 = np.append(new_order_top_100, np.where(cluster_labels==c)[0])

    #re-order
    all_metrics_arr_top_100_reordered = all_metrics_arr_top_100[:, new_order_top_100]
    top_region_inds_reordered = np.append(top_region_inds[0], top_region_inds[1:][new_order_top_100])
    score_arr_top_100 = score_arr_top_100[new_order_top_100]
    
    metric_sum_top_100 = np.sum(all_metrics_arr_top_100_reordered, axis=0)
    metrics_bar_scaled = all_metrics_arr_top_100_reordered / metric_sum_top_100 * score_arr_top_100 

    metrics_bar_scaled_inds = np.argsort(-metrics_bar_scaled, axis=0)
    metrics_bar_scaled_vals = -np.sort(-metrics_bar_scaled, axis=0)

    #all string original and destination genes
    all_arr_ori = all_arr_ori[:, new_order_top_100]
    all_arr_dest = all_arr_dest[:, new_order_top_100]

    #plotting
    fig, axes = plt.subplots(ncols=2, figsize=(16, 20))
    fig.patch.set_facecolor('white')
    plt.subplots_adjust(wspace=0)
    axes = axes.ravel()

    #sideways bar plot
    for i in range(11):
        axes[0].barh(np.flip(np.arange(100)), metrics_bar_scaled_vals[i], left=np.sum(metrics_bar_scaled_vals[:i], axis=0), 
                 color=colors[metrics_bar_scaled_inds[i]], height=1, align='center')

    #label top comps
    top_factors = []
    for i in range(len(metrics_bar_scaled_vals[0])):
        start_val = 0
        for metric, val in zip(metrics_bar_scaled_inds[:, i], metrics_bar_scaled_vals[:, i]):
            if ((metric==0) & (val>0)):
                DHS_i_rank = rankdata(-DHS_ooe[top_region_inds_reordered[i+1]])
                top_comp = np.argsort((DHS_roi_rank+DHS_i_rank)/2)[0]
                if val / len(DHS_component_order[top_comp]) > 0.012:
                    axes[0].text(start_val + val/10, 
                             len(metrics_bar_scaled_vals[0])-1.35-i, 
                             DHS_component_order[top_comp], 
                             fontsize=9)
                    top_factors.append(DHS_component_order[top_comp])
            elif ((metric==1) & (val>0)):
                TF_i_rank = rankdata(-TF_ooe[top_region_inds_reordered[i+1]])
                top_comp = np.argsort((TF_roi_rank+TF_i_rank)/2)[0]
                if val / len(TF_component_order[top_comp]) > 0.012:
                    axes[0].text(start_val + val/10, 
                             len(metrics_bar_scaled_vals[0])-1.35-i, 
                             TF_component_order[top_comp], 
                             fontsize=9)
                    top_factors.append(TF_component_order[top_comp])
            elif ((metric==2) & (val>0)):
                epilogos_i_mean_rank = rankdata(-epilogos_ooe[top_region_inds_reordered[i+1]][:18])
                epilogos_i_max_rank = rankdata(-epilogos_ooe[top_region_inds_reordered[i+1]][18:])
                top_mean_comp = np.argsort((epilogos_roi_mean_rank + epilogos_i_mean_rank)/2)[0]
                top_max_comp = np.argsort((epilogos_roi_max_rank + epilogos_i_max_rank)/2)[0]
                if val / len(epilogos_component_order[top_mean_comp] + '/' + epilogos_component_order[top_max_comp]) > 0.012:
                    axes[0].text(start_val + val/10, 
                             len(metrics_bar_scaled_vals[0])-1.35-i, 
                             epilogos_component_order[top_mean_comp] + '/' + epilogos_component_order[top_max_comp], 
                             fontsize=9)
                    top_factors.append(epilogos_component_order[top_mean_comp] + '/' + epilogos_component_order[top_max_comp])
            elif ((metric>2) & (val>0)):
                metric_ori_gene = geneID2name[all_arr_ori[metric-3][i]]
                metric_dest_gene = geneID2name[all_arr_dest[metric-3][i]]
                if val / len(metric_ori_gene + ' to ' + metric_dest_gene) > 0.0135:
                    axes[0].text(start_val + val/10, 
                             len(metrics_bar_scaled_vals[0])-1.35-i, 
                             metric_ori_gene + ' to ' + metric_dest_gene, 
                             fontsize=9)
                    top_factors.append(metric_ori_gene + ' to X')
            start_val+=val
        axes[1].text(0.01, 
                     (len(metrics_bar_scaled_vals[0])-0.75-i)/len(metrics_bar_scaled_vals[0]), 
                     np.round(score_arr_top_100[i], 4), 
                     fontsize=9)

    axes[1].text(0.01, 1.005, 'Overall Score', fontsize=9)

    #Seperating Clusters with black lines
    height = len(metrics_bar_scaled_vals[0])-0.5
    for c in np.unique(cluster_labels)[:-1]: #not inclding line at bottom
        num_in_cluster = len(cluster_labels[cluster_labels==c])
        height -= num_in_cluster
        axes[0].plot([0, 1], [height, height], c='r', linewidth=2)

    #look up table
    top_factors = np.unique(top_factors)
    factor_bar_arr = []
    for factor in top_factors:
        factor_l = factor.split('/')
        if len(factor_l) == 2: #epilogos
            if factor_l[0]==factor_l[1]:
                factor_bar_arr.append(epilogosfact2descr[factor_l[0]])
            else:
                factor_bar_arr.append(epilogosfact2descr[factor_l[0]] + ', with ' + epilogosfact2descr[factor_l[1]] + ' peaks')
        elif factor in DHSfact2descr: #DHS
            factor_bar_arr.append(DHSfact2descr[factor])
        elif factor in TFfact2descr: #TF motifs
            factor_bar_arr.append(factor + ': ' + TFfact2descr[factor])
        else:
            gene0, _ = factor.split(' to ')[0], factor.split(' to ')[1]
            factor_bar_arr.append('Gene ' + gene0 + ' Mapping to Gene X')
    factor_bar_arr = np.array(factor_bar_arr)
    lookup_table = np.hstack((top_factors.reshape((-1, 1)), factor_bar_arr.reshape((-1, 1))))
    t = axes[1].table(lookup_table, 
                  cellLoc='left',
                  colLabels=['Key', 'Descriptor'],
                  bbox=[0.2, 1-(len(factor_bar_arr)/20), 0.75, (len(factor_bar_arr)/20)], 
                  colWidths=[0.15, 0.35])
    t.auto_set_font_size(False)
    t.set_fontsize(9)

    axes[1].set_xticks([])
    axes[1].set_yticks([])

    axes[1].spines['right'].set_visible(False)
    # axes[1].spines['left'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['bottom'].set_visible(False)

    axes[0].set_yticks([])
    axes[0].set_ylim([-0.5, 99.5])
    axes[0].set_xticks([1])
    h = []

    for metric, c in zip(['DNase I Accessibility', 'TF Motifs', 'Chromatin States', 'Co-expression', 
                          'Co-evolved', 'Co-annotated', 'Interacts', 'Fused', 'Homologous', 
                          'Neighborhood', 'Co-cited'], colors):
        p = mpatches.Patch(color=c, label=metric)
        h.append(p)
    axes[1].legend(handles=h, loc='lower right', ncol=3)
    axes[0].set_xlabel('Overall Score')
    title = 'Scoring Breakdown for ROI %s (Mean Shift Clustering)' %(roi)
    axes[0].set_title(title)
    plt.savefig(figure_out_dir + title.replace(' ', '_') + '_CT' + today.strftime("%Y%m%d") + '.pdf', 
                bbox_inches='tight', 
                dpi=400, 
                facecolor="#FFFFFF", 
                edgecolor="#FFFFFF", 
                transparent=False)
    return

def make_reduced_barplot(all_metrics_arr_top_100, top_region_inds, all_arr_ori, colors, 
                         DHS_roi_rank, TF_roi_rank, epilogos_roi_mean_rank, epilogos_roi_max_rank, DHS_ooe, 
                         TF_ooe, epilogos_ooe, DHS_component_order, TF_component_order, epilogos_component_order, 
                         geneID2name, DHSfact2descr, TFfact2descr, epilogosfact2descr, roi, figure_out_dir):
    
    #sideways bar plots (average of each kmeans cluster)
    #cluster rows
    cluster = MeanShift()
    cluster.fit(all_metrics_arr_top_100.T)
    clusters_for_each_region = cluster.labels_

    # new_order_best_of_clusters = []
    cluster_percentages = np.array([])
    score_arr_top = np.array([])
    all_metrics_arr_top = []
    regions_per_cluster = [] #for labeling top factors
    for cluster_ind in np.unique(clusters_for_each_region):
        regions_in_cluster = np.where(clusters_for_each_region==cluster_ind)[0]
        regions_per_cluster.append(regions_in_cluster)
        cluster_percentages = np.append(cluster_percentages, int(len(regions_in_cluster)/len(clusters_for_each_region)*100)) #as percent
        cluster_scores = all_metrics_arr_top_100[:, regions_in_cluster]
        cluster_scores_mean = np.mean(cluster_scores, axis=1)
        score = 1 - np.prod(1 - cluster_scores_mean)
        score_arr_top = np.append(score_arr_top, score)
        all_metrics_arr_top.append(cluster_scores_mean)
    all_metrics_arr_top = np.array(all_metrics_arr_top).T
    regions_per_cluster = np.array(regions_per_cluster, dtype=object)

    #sorting
    # sort_order = np.argsort(-score_arr_top)
    sort_order = np.argsort(-cluster_percentages)
    all_metrics_arr_top = all_metrics_arr_top[:, sort_order]
    cluster_percentages = cluster_percentages[sort_order]
    score_arr_top = score_arr_top[sort_order]
    regions_per_cluster = regions_per_cluster[sort_order]

    metric_sum_top = np.sum(all_metrics_arr_top, axis=0)

    metrics_bar_scaled = all_metrics_arr_top / metric_sum_top * score_arr_top

    metrics_bar_scaled_inds = np.argsort(-metrics_bar_scaled, axis=0)
    metrics_bar_scaled_vals = -np.sort(-metrics_bar_scaled, axis=0)

    #plotting
    fig, axes = plt.subplots(ncols=2, figsize=(16, len(sort_order)*0.52))
    fig.patch.set_facecolor('white')
    plt.subplots_adjust(wspace=0.3)
    axes = axes.ravel()

    #sideways bar plot
    for i in range(11):
        axes[0].barh(np.flip(np.arange(len(score_arr_top))), 
                     metrics_bar_scaled_vals[i], 
                     left=np.sum(metrics_bar_scaled_vals[:i], axis=0), 
                     color=colors[metrics_bar_scaled_inds[i]], height=1, align='center')

    #label top comps
    top_factors = []
    for i in range(len(metrics_bar_scaled_vals[0])):
        start_val = 0
        regions_in_cluster = regions_per_cluster[i]
        for metric, val in zip(metrics_bar_scaled_inds[:, i], metrics_bar_scaled_vals[:, i]):
            if ((metric==0) & (val>0)):
                DHS_i_rank = rankdata(-np.mean(DHS_ooe[top_region_inds[1:][regions_in_cluster]], axis=0))
                top_comp = np.argsort((DHS_roi_rank+DHS_i_rank)/2)[0]
                if val / len(DHS_component_order[top_comp]) > 0.012:
                    axes[0].text(start_val + val/10, 
                             len(metrics_bar_scaled_vals[0])-1.35-i, 
                             DHS_component_order[top_comp], 
                             fontsize=9)
                    top_factors.append(DHS_component_order[top_comp])
            elif ((metric==1) & (val>0)):
                TF_i_rank = rankdata(-np.mean(TF_ooe[top_region_inds[1:][regions_in_cluster]], axis=0))
                top_comp = np.argsort((TF_roi_rank+TF_i_rank)/2)[0]
                if val / len(TF_component_order[top_comp]) > 0.012:
                    axes[0].text(start_val + val/10, 
                             len(metrics_bar_scaled_vals[0])-1.35-i, 
                             TF_component_order[top_comp], 
                             fontsize=9)
                    top_factors.append(TF_component_order[top_comp])
            elif ((metric==2) & (val>0)):
                epilogos_i_mean_rank = rankdata(-np.mean(epilogos_ooe[top_region_inds[1:][regions_in_cluster]][:, :18], axis=0))
                epilogos_i_max_rank = rankdata(-np.mean(epilogos_ooe[top_region_inds[1:][regions_in_cluster]][:, 18:], axis=0))
                top_mean_comp = np.argsort((epilogos_roi_mean_rank + epilogos_i_mean_rank)/2)[0]
                top_max_comp = np.argsort((epilogos_roi_max_rank + epilogos_i_max_rank)/2)[0]
                if val / len(epilogos_component_order[top_mean_comp] + '/' + epilogos_component_order[top_max_comp]) > 0.012:
                    axes[0].text(start_val + val/10, 
                             len(metrics_bar_scaled_vals[0])-1.35-i, 
                             epilogos_component_order[top_mean_comp] + '/' + epilogos_component_order[top_max_comp], 
                             fontsize=9)
                    top_factors.append(epilogos_component_order[top_mean_comp] + '/' + epilogos_component_order[top_max_comp])
            elif ((metric>2) & (val>0)):
                genes, counts = np.unique(all_arr_ori[:, regions_in_cluster][metric-3].astype(str), return_counts=True)
                genes_ordered = genes[np.argsort(-counts)]
                if genes_ordered[0] != 'None':
                    metric_ori_gene = geneID2name[genes_ordered[0]]
                else:
                    metric_ori_gene = geneID2name[genes_ordered[1]]
                if val / len(metric_ori_gene + ' to X') > 0.0135:
                    axes[0].text(start_val + val/10, 
                             len(metrics_bar_scaled_vals[0])-1.35-i, 
                             metric_ori_gene + ' to X', 
                             fontsize=9)
                    top_factors.append(metric_ori_gene + ' to X')
            start_val+=val

    #fusion scores and percent of top 100 in each cluster
    for i in range(len(metrics_bar_scaled_vals[0])):
        axes[0].text(1.01, 
                     len(metrics_bar_scaled_vals[0])-1.35-i, 
                     np.round(score_arr_top[i], 4), 
                     fontsize=9)
        axes[0].text(1.18, 
                     len(metrics_bar_scaled_vals[0])-1.35-i, 
                     str(int(cluster_percentages[i])) + '%', 
                     fontsize=9)

    axes[0].text(1.01, len(metrics_bar_scaled_vals[0]), 'Overall Score', fontsize=9)
    axes[0].text(1.18, len(metrics_bar_scaled_vals[0]), 'Cluster Percentage', fontsize=9)

    #look up table
    top_factors = np.unique(top_factors)
    factor_bar_arr = []
    for factor in top_factors:
        factor_l = factor.split('/')
        if len(factor_l) == 2: #epilogos
            if factor_l[0]==factor_l[1]:
                factor_bar_arr.append(epilogosfact2descr[factor_l[0]])
            else:
                factor_bar_arr.append(epilogosfact2descr[factor_l[0]] + ', with ' + epilogosfact2descr[factor_l[1]] + ' peaks')
        elif factor in DHSfact2descr: #DHS
            factor_bar_arr.append(DHSfact2descr[factor])
        elif factor in TFfact2descr: #TF motifs
            factor_bar_arr.append(factor + ': ' + TFfact2descr[factor])
        else:
            gene0, gene1 = factor.split(' to ')[0], factor.split(' to ')[1]
            factor_bar_arr.append('Gene ' + gene0 + ' Mapping to Gene X')
    factor_bar_arr = np.array(factor_bar_arr)
    lookup_table = np.hstack((top_factors.reshape((-1, 1)), factor_bar_arr.reshape((-1, 1))))
    t = axes[1].table(lookup_table, 
                      cellLoc='left',
                      colLabels=['Key', 'Descriptor'],
                      bbox=[0.15, 0.3, 0.85, 0.7], 
                      colWidths=[0.15, 0.35])
    t.auto_set_font_size(False)
    t.set_fontsize(9)

    axes[1].set_xticks([])
    axes[1].set_yticks([])

    axes[1].spines['right'].set_visible(False)
    axes[1].spines['left'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['bottom'].set_visible(False)

    axes[0].set_yticks([])
    axes[0].set_ylim([-7, len(all_metrics_arr_top[0])-0.5])
    axes[0].set_xticks([1])
    h = []

    for metric, c in zip(['DNase I Accessibility', 'TF Motifs', 'Chromatin States', 'Co-expression', 
                          'Co-evolved', 'Co-annotated', 'Interacts', 'Fused', 'Homologous', 
                          'Neighborhood', 'Co-cited'], colors):
        p = mpatches.Patch(color=c, label=metric)
        h.append(p)
    axes[0].legend(handles=h, loc='lower center', ncol=3)
    axes[0].set_xlabel('Overall Score')
    title = 'Scoring Breakdown for ROI %s (Average, Mean Shift Cluster)' %(roi)
    axes[0].set_title(title)
    plt.savefig(figure_out_dir + title.replace(' ', '_') + '_CT' + today.strftime("%Y%m%d") + '.pdf', 
                bbox_inches='tight', 
                dpi=400, 
                facecolor="#FFFFFF", 
                edgecolor="#FFFFFF", 
                transparent=False)
    return


def main(gz_file_path, genomic_location, DHS_hdf5_path, TF_hdf5_path, epilogos_hdf5_path, coex_hdf5_path, cooccur_hdf5_path, database_hdf5_path, exper_hdf5_path,
         fusion_hdf5_path, homology_hdf5_path, neigh_hdf5_path, textmin_hdf5_path, regions_path, DHS_ooe_path, TF_ooe_path, epilogos_ooe_path,
         DHS_mapped_path, TF_mapped_path, epilogos_mapped_path, DHSfact2descr_path,  TFfact2descr_path, epifact2descr_path, geneID2name_path, out_bed_path, 
         figure_out_dir):

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

    for r in range(len(roi)):
        make_barplot(top_region_inds[r], metrics_bar_scaled_vals[:, r], metrics_bar_scaled_inds[:, r], 
                     DHS_roi_rank[r], TF_roi_rank[r], epilogos_roi_mean_rank[r], epilogos_roi_max_rank[r], 
                     DHS_ooe, TF_ooe, epilogos_ooe, DHS_component_order, TF_component_order, 
                     epilogos_component_order, all_arr_ori[:, r], all_arr_dest[:, r], geneID2name, 
                     score_arr_top_100[r], DHSfact2descr, TFfact2descr, epilogosfact2descr, roi[r], colors, 
                     figure_out_dir)

        make_cluster_barplot(all_metrics_arr_top_100[:, r], top_region_inds[r], score_arr_top_100[r], all_arr_ori[:, r],
                             all_arr_dest[:, r], colors, DHS_roi_rank[r], TF_roi_rank[r], epilogos_roi_mean_rank[r], 
                             epilogos_roi_max_rank[r], DHS_ooe, TF_ooe, epilogos_ooe, DHS_component_order, 
                             TF_component_order, epilogos_component_order, geneID2name, DHSfact2descr, TFfact2descr, 
                             epilogosfact2descr, roi[r], figure_out_dir)

        make_reduced_barplot(all_metrics_arr_top_100[:, r], top_region_inds[r], all_arr_ori[:, r],
                             colors, DHS_roi_rank[r], TF_roi_rank[r], epilogos_roi_mean_rank[r], 
                             epilogos_roi_max_rank[r], DHS_ooe, TF_ooe, epilogos_ooe, DHS_component_order, 
                             TF_component_order, epilogos_component_order, geneID2name, DHSfact2descr, TFfact2descr, 
                             epilogosfact2descr, roi[r], figure_out_dir)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12],
         sys.argv[13], sys.argv[14], sys.argv[15], sys.argv[16], sys.argv[17], sys.argv[18], sys.argv[19], sys.argv[20], sys.argv[21], sys.argv[22], sys.argv[23], 
         sys.argv[24], sys.argv[25], sys.argv[26])
