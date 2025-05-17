import numpy as np
from functools import partial
from multiprocessing import Pool
from contextlib import closing
import sys
import h5py
from collections import Counter
import re
from region_comparison_combine_scores_CT20220711.py import combine_scores

def curate_descriptor_terms(genes_bed_path):
    #curate a list of descriptor terms
    f = open(genes_bed_path, 'r')
    descriptor_terms = np.array([])
    for line in f.readlines():
        line = line.strip().split('\t')
        descriptor = line[4].lower()
        descriptor = descriptor.replace('(', '').replace(')', '')
        descriptor = np.array(descriptor.split())
        for i in range(len(descriptor)):
            descriptor[i] = descriptor[i].strip(',')
        descriptor_terms = np.append(descriptor_terms, descriptor)
    descriptor_terms = np.unique(descriptor_terms)
    f.close()

    stop_words = ['of', 'the', 'a', 'protein', 'and', 'like', 'with', 'domain', 'containing', 'member', 
                  'associated', 'family', 'non', 'paired', 'one', 'two', 'three', 'four', 'five', 'six', 
                  'seven', 'eight', 'nine', 'ten', 'class', 'related', 'bound', 'type', 'pro', '-']
    unwanted_chars = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '-', ',', '/', '.']
    for descriptor in descriptor_terms:
        char_count_dict = Counter(descriptor)
        unwanted_sum = 0
        for char in unwanted_chars:
            unwanted_sum += char_count_dict[char]
        if ((descriptor.isnumeric()) | 
            (descriptor in stop_words) | 
            (unwanted_sum >= len(descriptor) / 2) | 
            (len(descriptor) < 2)):
            descriptor_terms = descriptor_terms[descriptor_terms!=descriptor]
    return descriptor_terms

def count_terms(descriptor_terms, DHS_scaled_vals, DHS_similar_regions, TF_scaled_vals, TF_similar_regions, epilogos_scaled_vals, epilogos_similar_regions, coex_scaled_vals, 
		coex_similar_regions, cooccur_scaled_vals, cooccur_similar_regions, database_scaled_vals, database_similar_regions, exper_scaled_vals, exper_similar_regions,
		fusion_scaled_vals, fusion_similar_regions, homology_scaled_vals, homology_similar_regions, neigh_scaled_vals, neigh_similar_regions, textmin_scaled_vals, 
		textmin_similar_regions, mapped_genes, roi):

    term_roi_counts = np.zeros(len(descriptor_terms), dtype=int)

    score_arr, _ = combine_scores(DHS_scaled_vals[roi], DHS_similar_regions[roi], TF_scaled_vals[roi], TF_similar_regions[roi], epilogos_scaled_vals[roi],
                  		  epilogos_similar_regions[roi], coex_scaled_vals[roi], coex_similar_regions[roi], cooccur_scaled_vals[roi], cooccur_similar_regions[roi],
              		          database_scaled_vals[roi], database_similar_regions[roi], exper_scaled_vals[roi], exper_similar_regions[roi], fusion_scaled_vals[roi], fusion_similar_regions[roi],
	                          homology_scaled_vals[roi], homology_similar_regions[roi], neigh_scaled_vals[roi], neigh_similar_regions[roi], textmin_scaled_vals[roi], textmin_similar_regions[roi])

    #finding top inds
    rec_inds = np.argsort(-score_arr)[:100]
    top_region_inds = np.append(roi, rec_inds)

    #find descriptors for mapped genes
    roi_descriptions = ''
    for i in top_region_inds:
        line = re.split(';|\t|\|', mapped_genes[i].strip())
        descriptions = line[7::5]
        roi_descriptions += ''.join(descriptions) + ' '

    #processing gene descriptions
    roi_descriptions = roi_descriptions.replace('(', '').replace(')', '').lower().split()
    for descriptor in roi_descriptions:
        descriptor = descriptor.strip(',')
        ind_l = np.where(descriptor_terms==descriptor)[0]
        if len(ind_l) == 1:
            term_roi_counts[ind_l[0]] += 1

    return term_roi_counts

def ooe_calc(term_roi_counts):
    #calculate ooe for each term in each roi
    totals_per_term = np.sum(term_roi_counts, axis=0)
    total_count = np.sum(totals_per_term)
    expected = totals_per_term / total_count

    totals_per_roi = np.sum(term_roi_counts, axis=1)
    observed = (term_roi_counts.T / totals_per_roi).T

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
         fusion_hdf5_path, homology_hdf5_path, neigh_hdf5_path, textmin_hdf5_path, genes_bed_path, mapped_genes_path, term_ooe_path, 
         terms_path, num_cpus):
    #load the data
    print('Loading data...', flush=True)
    ###Load rank and ind arrays
    #DHS counts
    DHS_scaled_vals, DHS_similar_regions = load_hdf5_files(DHS_hdf5_path, roi)
    #TF motif counts
    TF_scaled_vals, TF_similar_regions = load_hdf5_files(TF_hdf5_path, roi)
    #epilogos score sum and max
    epilogos_scaled_vals, epilogos_similar_regions = load_hdf5_files(epilogos_hdf5_path, roi)
    #coexpression STRING scores
    coex_scaled_vals, coex_similar_regions = load_hdf5_files(coex_hdf5_path, roi)
    #cooccurence STRING scores
    cooccur_scaled_vals, cooccur_similar_regions = load_hdf5_files(cooccur_hdf5_path, roi)
    #database STRING scores
    database_scaled_vals, database_similar_regions = load_hdf5_files(database_hdf5_path, roi)
    #experiments STRING scores
    exper_scaled_vals, exper_similar_regions = load_hdf5_files(exper_hdf5_path, roi)
    #fusion STRING scores
    fusion_scaled_vals, fusion_similar_regions = load_hdf5_files(fusion_hdf5_path, roi)
    #homology STRING scores
    homology_scaled_vals, homology_similar_regions = load_hdf5_files(homology_hdf5_path, roi)
    #neighborhood STRING scores
    neigh_scaled_vals, neigh_similar_regions = load_hdf5_files(neigh_hdf5_path, roi)
    #textmining STRING scores
    textmin_scaled_vals, textmin_similar_regions = load_hdf5_files(textmin_hdf5_path, roi)


    DHS_hdf5.close()
    TF_hdf5.close()
    epilogos_hdf5.close()

    f = open(mapped_genes_path, 'r')
    mapped_genes = f.readlines()
    f.close()

    #generate list of descriptor terms
    print('Generating term list...', flush=True)
    descriptor_terms = curate_descriptor_terms(genes_bed_path)

    #split on cluster and find similar regions based on cosine similarity
    print('Counting terms for each region...', flush=True)
    func = partial(count_terms, descriptor_terms, DHS_scaled_vals, DHS_similar_regions, TF_scaled_vals, TF_similar_regions, epilogos_scaled_vals, epilogos_similar_regions, coex_scaled_vals, 
                   coex_similar_regions, cooccur_scaled_vals, cooccur_similar_regions, database_scaled_vals, database_similar_regions, exper_scaled_vals, exper_similar_regions,
                   fusion_scaled_vals, fusion_similar_regions, homology_scaled_vals, homology_similar_regions, neigh_scaled_vals, neigh_similar_regions, textmin_scaled_vals, 
                   textmin_similar_regions, mapped_genes)
    num_regions = len(DHS_scaled_vals)
    with closing(Pool(num_cpus)) as pool:
        term_roi_counts = pool.map(func, np.arange(num_regions))
    pool.join()
    term_roi_counts = np.array(term_roi_counts)

    #ooe of each term
    ooe = ooe_calc(term_roi_counts)

    #save the terms and the term counts for each region
    ooe_f = h5py.File(term_ooe_path, 'w')
    dset = ooe_f.create_dataset('term_roi_counts', data=term_roi_counts)
    dset = ooe_f.create_dataset('term_roi_ooe', data=ooe)
    ooe_f.close()
    
    np.save(terms_path, descriptor_terms)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9] sys.argv[10] sys.argv[11] sys.argv[12] sys.argv[13] sys.argv[14] sys.argv[15] int(sys.argv[16]))
