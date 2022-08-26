import numpy as np
import sys
import h5py
from sklearn.metrics.pairwise import cosine_similarity
from scipy.cluster import hierarchy

def main(genesetooe_sig_hdf5_path, out_path):
    #load data
    geneset_hdf5 = h5py.File(genesetooe_sig_hdf5_path, 'r')
    gene_set_sig_log2ooe = geneset_hdf5['sig_gene_set_roi_log2ooe'][:]
    gene_set_sig_counts = geneset_hdf5['sig_gene_set_roi_counts'][:]
    geneset_hdf5.close()

    #create distance matrix
    gene_set_cosine_sim_m = cosine_similarity(gene_set_sig_log2ooe.T)
    #linkage
    average_cluster = hierarchy.average(gene_set_cosine_sim_m)
    #get the optimal ordering
    average_cluster_reordered = hierarchy.optimal_leaf_ordering(average_cluster, gene_set_cosine_sim_m)
    #get the terms in their clustering order
    gene_set_order = hierarchy.leaves_list(average_cluster_reordered)

    #save the term order
    np.save(out_path, gene_set_order)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
