import numpy as np
import sys
import os
import time
from datetime import date
import h5py
import tabix
from scipy.stats import rankdata
import pickle
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import cm, colorbar
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
today = date.today()

def make_spectrum_plot(roi, gene_set_order, inds_clustered, sig_file_path, all_file_path, 
                       regions_path, figure_out_dir):
    #gene set spectrum plot
    hdf5_file = h5py.File(sig_file_path, 'r')
    gene_set_sig_log2ooe = hdf5_file['sig_gene_set_roi_log2ooe'][roi][inds_clustered][gene_set_order]
    gene_set_sig_counts = hdf5_file['sig_gene_set_roi_counts'][roi][inds_clustered][gene_set_order]
    hdf5_file.close()
    
    hdf5_file = h5py.File(all_file_path, 'r')
    gene_sets = np.char.decode(hdf5_file['gene_sets'][:], 'utf-8')[inds_clustered][gene_set_order]
    hdf5_file.close()
    
    colors = cm.get_cmap('nipy_spectral', len(gene_set_order))
    
    #max over a given window size
    window_size = 20
    gene_set_sig_log2ooe_windowed = np.array([])
    for i in range(len(gene_set_sig_log2ooe)):
        gene_set_sig_log2ooe_windowed = np.append(gene_set_sig_log2ooe_windowed, 
                                                  np.max(gene_set_sig_log2ooe[max(i-int(window_size/2), 0):
                                                                              min(i+int(window_size/2), 
                                                                                  len(gene_set_sig_log2ooe))]))
    #smoothing with mean
    gene_set_sig_log2ooe_windowed_smoothed = np.array([])
    for i in range(len(gene_set_sig_log2ooe_windowed)):
        gene_set_sig_log2ooe_windowed_smoothed = np.append(gene_set_sig_log2ooe_windowed_smoothed, 
                                                           np.mean(gene_set_sig_log2ooe_windowed[max(i-int(window_size/2), 0):
                                                                                                 min(i+int(window_size/2), 
                                                                                                     len(gene_set_sig_log2ooe_windowed))]))
    
    cluster_inds_sorted = np.where(gene_set_sig_log2ooe!=0)[0]
    gene_sets_labels_ordered = gene_sets[cluster_inds_sorted]
    gene_sets_vals_ordered = gene_set_sig_log2ooe[cluster_inds_sorted]
    
    gene_sets_to_label = np.array([])
    cluster_inds_to_label = np.array([])
    vals_to_label = np.array([])
    for cluster_ind, label, val in zip(cluster_inds_sorted, gene_sets_labels_ordered, gene_sets_vals_ordered):
        label = ' '.join(label.split('_'))
        if len(vals_to_label) == 0:
            gene_sets_to_label = np.append(gene_sets_to_label, label)
            cluster_inds_to_label = np.append(cluster_inds_to_label, cluster_ind)
            vals_to_label = np.append(vals_to_label, val)
        elif ((cluster_ind - cluster_inds_to_label[-1] <= 100) & (vals_to_label[-1] < val)):
            gene_sets_to_label[-1] = label
            cluster_inds_to_label[-1] = cluster_ind
            vals_to_label[-1] = val
        elif cluster_ind - cluster_inds_to_label[-1] > 100:
            gene_sets_to_label = np.append(gene_sets_to_label, label)
            cluster_inds_to_label = np.append(cluster_inds_to_label, cluster_ind)
            vals_to_label = np.append(vals_to_label, val)
            
    #trimming terms
    for i in range(len(gene_sets_to_label)):
        gene_sets_to_label[i] = gene_sets_to_label[i].replace('WP ', '').replace('HALLMARK ', '').replace('KEGG ', '').replace('REACTOME ', '').replace('BIOCARTA ', '')

    #plot
    fig, axes = plt.subplots(nrows=2, gridspec_kw={'height_ratios': [4, 1]}, figsize=(12, 2))
    fig.patch.set_facecolor('white')
    axes = axes.ravel()
    plt.subplots_adjust(hspace=0)

    axes[0].fill_between(range(len(gene_set_sig_log2ooe_windowed_smoothed)), 
                         gene_set_sig_log2ooe_windowed_smoothed, color='k')

    colorbar.ColorbarBase(axes[1], cmap=colors, orientation='horizontal').set_ticks(cluster_inds_to_label/len(gene_sets))
    axes[0].set_xlim([0, len(gene_sets)])
    axes[0].set_xticks([])
    axes[0].set_ylabel('Windowed Enrichment')
    axes[0].set_ylim([0, np.max(gene_set_sig_log2ooe_windowed)*1.1])
    axes[1].set_xticklabels(gene_sets_to_label, rotation=45, ha='right')

    for xlabel, size in zip(axes[1].xaxis.get_ticklabels(), vals_to_label):
        xlabel.set_fontsize(size*2.5)


    f = open(regions_path, 'r')
    regions = f.readlines()
    f.close()
    regions = regions[roi].strip().split('\t')
    genome_coords = regions[0] + ':' + regions[1] + '-' + regions[2]

    title = 'Smoothed Windowed Enrichment of Significant OOE Gene Sets for ROI %s\n(%s) and Top 100 Regions' %(roi, genome_coords)
    fig.suptitle(title, y=1.08)
    plt.savefig(figure_out_dir + title.replace(' ', '_').replace(':', '').replace('(', '').replace(')', '').replace('\n', '_') + '_CT' + today.strftime("%Y%m%d") + '.pdf', 
                bbox_inches='tight', 
                dpi=400, 
                facecolor="#FFFFFF", 
                edgecolor="#FFFFFF", 
                transparent=False)
    return

def main(gz_file_path, genomic_location, regions_path, gene_set_order_path, sig_file_path, all_file_path, figure_out_dir):

    #find regions of interest    
    gz_tb = tabix.open(gz_file_path)
    roi = np.array(list(gz_tb.querys(genomic_location)))[:, 3].astype(int)
    print('roi = ' + str(roi), flush=True)

    #gene set plots
    gene_set_order = np.load(gene_set_order_path)['gene_set_order']
    inds_clustered = np.load(gene_set_order_path)['inds_clustered']

    for r in range(len(roi)):
        make_spectrum_plot(roi[r], gene_set_order, inds_clustered, sig_file_path, 
                           all_file_path, regions_path, figure_out_dir)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
