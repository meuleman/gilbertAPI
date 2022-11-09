import numpy as np
import os
import sys
import h5py
from scipy.stats import rankdata

def load_STRING_hdf5_files_multiple(hdf5_path):
    hdf5_file = h5py.File(hdf5_path, 'r')
    scaled_vals = hdf5_file['top_1000_scaled_vals'][:]
    hdf5_file.close()
    return scaled_vals

def rank_regions(arr):
    #rank each region
    ranked_regions = rankdata(arr, axis=0)

    #find the number of regions with 0 values for each component
    comp_num_zero = len(arr) - np.count_nonzero(arr, axis=0)

    #adjusting/resetting ranks
    ranked_regions = ranked_regions - comp_num_zero
    ranked_regions = np.where(ranked_regions < 0, 0, ranked_regions)

    return ranked_regions

def main(coex_hdf5_path, cooccur_hdf5_path, database_hdf5_path, exper_hdf5_path, fusion_hdf5_path, homology_hdf5_path, neigh_hdf5_path, 
         textmin_hdf5_path, out_path):
    #load in individual STRING scores for all regions
    #coexpression STRING scores
    coex_scaled_vals = load_STRING_hdf5_files_multiple(coex_hdf5_path)
    #cooccurence STRING scores
    cooccur_scaled_vals = load_STRING_hdf5_files_multiple(cooccur_hdf5_path)
    #database STRING scores
    database_scaled_vals = load_STRING_hdf5_files_multiple(database_hdf5_path)
    #experiments STRING scores
    exper_scaled_vals = load_STRING_hdf5_files_multiple(exper_hdf5_path)
    #fusion STRING scores
    fusion_scaled_vals = load_STRING_hdf5_files_multiple(fusion_hdf5_path)
    #homology STRING scores
    homology_scaled_vals = load_STRING_hdf5_files_multiple(homology_hdf5_path)
    #neighborhood STRING scores
    neigh_scaled_vals = load_STRING_hdf5_files_multiple(neigh_hdf5_path)
    #textmining STRING scores
    textmin_scaled_vals = load_STRING_hdf5_files_multiple(textmin_hdf5_path)

    #averaging top 1000 individual STRING metric scores for each region
    regional_average_scores = []
    for scaled_vals in [coex_scaled_vals, cooccur_scaled_vals, database_scaled_vals, exper_scaled_vals, 
                        fusion_scaled_vals, homology_scaled_vals, neigh_scaled_vals, textmin_scaled_vals]:
        regional_average_scores.append(np.mean(scaled_vals, axis=1))
    regional_average_scores = np.array(regional_average_scores).T

    #ranking regions
    ranked_regions = rank_regions(regional_average_scores)

    #max ranks for each individual STRING metric
    max_ranks = np.max(ranked_regions, axis=0)

    #ordered metrics
    STRING_comps = ['Co-expression', 'Co-occurrence', 'Database', 'Experiment', 'Fusion', 'Homology', 
                    'Neighborhood', 'Textmining']

    #save ooe array
    out_f = h5py.File(out_path, 'w')
    dset = out_f.create_dataset('ranked_regions', data=ranked_regions)
    dset = out_f.create_dataset('max_ranks', data=max_ranks)
    dset = out_f.create_dataset('comp_order', data=np.char.encode(STRING_comps, 'utf-8'))
    out_f.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9])
