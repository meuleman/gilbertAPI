import numpy as np
import sys
import h5py
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist

def main(genesetooe_sig_hdf5_path, out_path):
    #load data
    print('Loading data...', flush=True)
    geneset_hdf5 = h5py.File(genesetooe_sig_hdf5_path, 'r')
    gene_set_sig_log2ooe = geneset_hdf5['sig_gene_set_roi_log2ooe'][:]
    geneset_hdf5.close()

    #just sets that are significant in at least 1 region
    sum_across_regions = np.sum(gene_set_sig_log2ooe, axis=0)
    inds_to_cluster = np.where(sum_across_regions>0)[0]
    gene_set_sig_log2ooe = gene_set_sig_log2ooe[:, inds_to_cluster]
    
    print('Distance matrix...', flush=True)
    #create distance matrix
    sigregionXgs_cospdist_m = pdist(gene_set_sig_log2ooe.T, metric='cosine')
    print('Linkage...', flush=True)
    #linkage
    sigregionXgs_cospdist_averagelinkage_opt = hierarchy.linkage(sigregionXgs_cospdist_m, method='average', optimal_ordering=True)
    #get the terms in their clustering order
    gene_set_order = hierarchy.leaves_list(sigregionXgs_cospdist_averagelinkage_opt)

    print('Saving...', flush=True)
    #save the term order
    np.savez(out_path, gene_set_order=gene_set_order, inds_clustered=inds_to_cluster, Z=sigregionXgs_cospdist_averagelinkage_opt, dist=sigregionXgs_cospdist_m)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
