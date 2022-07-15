import numpy as np
import sys
import h5py
from sklearn.metrics.pairwise import cosine_similarity
from scipy.cluster import hierarchy

def main(descriptorooe_sig_hdf5_path, out_path):
    #load data
    descriptor_hdf5 = h5py.File(descriptorooe_sig_hdf5_path, 'r')
    term_sig_log2ooe = descriptor_hdf5['sig_term_roi_log2ooe'][:]
    term_sig_counts = descriptor_hdf5['sig_term_roi_counts'][:]
    descriptor_hdf5.close()

    #create distance matrix
    term_cosine_sim_m = cosine_similarity(term_sig_log2ooe.T)
    #linkage
    ward_cluster = hierarchy.ward(term_cosine_sim_m)
    #get the optimal ordering
    ward_cluster_reordered = hierarchy.optimal_leaf_ordering(ward_cluster, term_cosine_sim_m)
    #get the terms in their clustering order
    term_order = hierarchy.leaves_list(ward_cluster_reordered)

    #save the term order
    np.save(out_path, term_order)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
