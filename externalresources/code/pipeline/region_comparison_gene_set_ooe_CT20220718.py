import numpy as np
from functools import partial
from multiprocessing import Pool
from contextlib import closing
import sys
import h5py
import re
import pickle
from region_comparison_combine_scores_CT20220713 import combine_scores

def count_terms(gene_sets, gene2geneset, DHS_scaled_vals, DHS_similar_regions, TF_scaled_vals, TF_similar_regions, epilogos_scaled_vals, epilogos_similar_regions, coex_scaled_vals,
                coex_similar_regions, cooccur_scaled_vals, cooccur_similar_regions, database_scaled_vals, database_similar_regions, exper_scaled_vals, exper_similar_regions,
                fusion_scaled_vals, fusion_similar_regions, homology_scaled_vals, homology_similar_regions, neigh_scaled_vals, neigh_similar_regions, textmin_scaled_vals,
                textmin_similar_regions, mapped_genes, roi):

    gene_set_roi_counts = np.zeros(len(gene_sets), dtype=int)

    score_arr, _, __ = combine_scores(DHS_scaled_vals[roi], DHS_similar_regions[roi], TF_scaled_vals[roi], TF_similar_regions[roi], epilogos_scaled_vals[roi],
                                  epilogos_similar_regions[roi], coex_scaled_vals[roi], coex_similar_regions[roi], cooccur_scaled_vals[roi], cooccur_similar_regions[roi],
                                  database_scaled_vals[roi], database_similar_regions[roi], exper_scaled_vals[roi], exper_similar_regions[roi], fusion_scaled_vals[roi], fusion_similar_regions[roi],
                                  homology_scaled_vals[roi], homology_similar_regions[roi], neigh_scaled_vals[roi], neigh_similar_regions[roi], textmin_scaled_vals[roi], textmin_similar_regions[roi])

    #finding top inds
    rec_inds = np.argsort(-score_arr)[:100]
    top_region_inds = np.append(roi, rec_inds)

    #find gene names for mapped genes
    roi_genes = np.array([])
    for i in top_region_inds:
        line = re.split(';|\t|\|', mapped_genes[i].strip())
        genes = line[6::5]
        roi_genes = np.append(roi_genes, genes)
    
    #collecting gene sets
    roi_gene_sets = np.array([])
    for gene in roi_genes:
        if gene in gene2geneset:
            roi_gene_sets = np.append(roi_gene_sets, gene2geneset[gene])
    
    #counting gene sets
    for gene_set in roi_gene_sets:
        ind = np.where(gene_sets==gene_set)[0][0]
        gene_set_roi_counts[ind] += 1

    return gene_set_roi_counts

def ooe_calc(gene_set_roi_counts):
    #calculate ooe for each gene set in each roi
    totals_per_gene_set = np.sum(gene_set_roi_counts, axis=0)
    total_count = np.sum(totals_per_gene_set)
    expected = totals_per_gene_set / total_count

    totals_per_roi = np.sum(gene_set_roi_counts, axis=1)
    observed = (gene_set_roi_counts.T / totals_per_roi).T

    ooe = observed / expected
    return ooe

def load_hdf5_files(hdf5_path):
    #take in a file path and load scaled scores and similar regions
    hdf5_file = h5py.File(hdf5_path, 'r')
    scaled_vals = hdf5_file['top_1000_scaled_vals'][:]
    similar_regions = hdf5_file['top_1000_similar_regions'][:]
    hdf5_file.close()
    return scaled_vals, similar_regions

def main(DHS_hdf5_path, TF_hdf5_path, epilogos_hdf5_path, coex_hdf5_path, cooccur_hdf5_path, database_hdf5_path, exper_hdf5_path,
         fusion_hdf5_path, homology_hdf5_path, neigh_hdf5_path, textmin_hdf5_path, mapped_genes_path, gene2geneset_path, gene_set_ooe_path, 
         num_cpus):
    #load the data
    print('Loading data...', flush=True)
    ###Load rank and ind arrays
    #DHS counts
    DHS_scaled_vals, DHS_similar_regions = load_hdf5_files(DHS_hdf5_path)
    #TF motif counts
    TF_scaled_vals, TF_similar_regions = load_hdf5_files(TF_hdf5_path)
    #epilogos score sum and max
    epilogos_scaled_vals, epilogos_similar_regions = load_hdf5_files(epilogos_hdf5_path)
    #coexpression STRING scores
    coex_scaled_vals, coex_similar_regions = load_hdf5_files(coex_hdf5_path)
    #cooccurence STRING scores
    cooccur_scaled_vals, cooccur_similar_regions = load_hdf5_files(cooccur_hdf5_path)
    #database STRING scores
    database_scaled_vals, database_similar_regions = load_hdf5_files(database_hdf5_path)
    #experiments STRING scores
    exper_scaled_vals, exper_similar_regions = load_hdf5_files(exper_hdf5_path)
    #fusion STRING scores
    fusion_scaled_vals, fusion_similar_regions = load_hdf5_files(fusion_hdf5_path)
    #homology STRING scores
    homology_scaled_vals, homology_similar_regions = load_hdf5_files(homology_hdf5_path)
    #neighborhood STRING scores
    neigh_scaled_vals, neigh_similar_regions = load_hdf5_files(neigh_hdf5_path)
    #textmining STRING scores
    textmin_scaled_vals, textmin_similar_regions = load_hdf5_files(textmin_hdf5_path)

    #getting genes mapped to each region
    f = open(mapped_genes_path, 'r')
    mapped_genes = f.readlines()
    f.close()

    #gene2geneset dictionary
    pickle_in = open(gene2geneset_path, 'rb')
    gene2geneset = pickle.load(pickle_in)

    #finding the different gene sets
    gene_sets = np.array([])
    for gene in gene2geneset.keys():
        gene_sets = np.append(gene_sets, gene2geneset[gene])
    gene_sets = np.unique(gene_sets)

    #split on cluster and find similar regions based on cosine similarity
    print('Counting terms for each region...', flush=True)
    func = partial(count_terms, gene_sets, gene2geneset, DHS_scaled_vals, DHS_similar_regions, TF_scaled_vals, TF_similar_regions, epilogos_scaled_vals, epilogos_similar_regions, coex_scaled_vals,
                   coex_similar_regions, cooccur_scaled_vals, cooccur_similar_regions, database_scaled_vals, database_similar_regions, exper_scaled_vals, exper_similar_regions,
                   fusion_scaled_vals, fusion_similar_regions, homology_scaled_vals, homology_similar_regions, neigh_scaled_vals, neigh_similar_regions, textmin_scaled_vals,
                   textmin_similar_regions, mapped_genes)
    num_regions = len(DHS_scaled_vals)
    with closing(Pool(num_cpus)) as pool:
        gene_set_roi_counts = pool.map(func, np.arange(num_regions))
    pool.join()
    gene_set_roi_counts = np.array(gene_set_roi_counts)

    #ooe of each term
    ooe = ooe_calc(gene_set_roi_counts)

    #save the terms and the term counts for each region
    ooe_f = h5py.File(gene_set_ooe_path, 'w')
    dset = ooe_f.create_dataset('gene_set_roi_counts', data=gene_set_roi_counts)
    dset = ooe_f.create_dataset('gene_set_roi_ooe', data=ooe)
    dset = ooe_f.create_dataset('gene_sets', data=np.char.encode(gene_sets, 'utf-8'))
    ooe_f.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13], sys.argv[14], int(sys.argv[15]))
