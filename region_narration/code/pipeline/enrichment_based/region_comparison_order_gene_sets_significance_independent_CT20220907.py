import numpy as np
import sys
import h5py
from sklearn.metrics.pairwise import cosine_similarity
from scipy.cluster import hierarchy

def main(genesetooe_hdf5_path, out_path):
    #load data
    print('Loading data...', flush=True)
    geneset_hdf5 = h5py.File(genesetooe_hdf5_path, 'r')
    gene_set_log2ooe = geneset_hdf5['gene_set_roi_ooe'][:]
    geneset_hdf5.close()

    #just sets that are present in at least 1 region
    sum_across_regions = np.sum(gene_set_log2ooe, axis=0)
    inds_to_cluster = np.where(sum_across_regions>0)[0]
    gene_set_log2ooe = gene_set_log2ooe[:, inds_to_cluster]
    
    print('Distance matrix...', flush=True)
    #create distance matrix
    gene_set_cosine_sim_m = cosine_similarity(gene_set_log2ooe.T)
    print('Linkage...', flush=True)
    #linkage
    average_cluster = hierarchy.average(gene_set_cosine_sim_m)
    print('Optimal order...', flush=True)
    #get the optimal ordering
    average_cluster_reordered = hierarchy.optimal_leaf_ordering(average_cluster, gene_set_cosine_sim_m)
    #get the terms in their clustering order
    gene_set_order = hierarchy.leaves_list(average_cluster_reordered)

    print('Saving...', flush=True)
    #save the term order
    np.savez(out_path, gene_set_order=gene_set_order, inds_clustered=inds_to_cluster)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
