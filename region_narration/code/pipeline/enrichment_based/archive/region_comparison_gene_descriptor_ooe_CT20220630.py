import numpy as np
from functools import partial
from multiprocessing import Pool
from contextlib import closing
import sys
import h5py
from collections import Counter
import re

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

def count_terms(descriptor_terms, DHS_scaled_vals, DHS_similar_regions, TF_scaled_vals, TF_similar_regions, epilogos_scaled_vals, epilogos_similar_regions, mapped_genes, roi):
    term_roi_counts = np.zeros(len(descriptor_terms), dtype=int)
    #DHS scores
    DHS_arr = np.zeros(len(DHS_scaled_vals))
    DHS_arr[DHS_similar_regions[roi]] = DHS_scaled_vals[roi]
    #epilogos scores
    epilogos_arr = np.zeros(len(DHS_scaled_vals))
    epilogos_arr[epilogos_similar_regions[roi]] = epilogos_scaled_vals[roi]
    #TF motif scores
    TF_arr = np.zeros(len(DHS_scaled_vals))
    TF_arr[TF_similar_regions[roi]] = TF_scaled_vals[roi]

    #scoring all metrics
    score_arr = 1 - ((1 - epilogos_arr) * (1 - DHS_arr) * (1 - TF_arr))

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

def main(DHS_hdf5_path, TF_hdf5_path, epilogos_hdf5_path, genes_bed_path, mapped_genes_path, term_ooe_path, terms_path, num_cpus):
    #load the data
    print('Loading data...', flush=True)
    DHS_hdf5 = h5py.File(DHS_hdf5_path, 'r')
    DHS_scaled_vals = DHS_hdf5['top_1000_scaled_vals'][:]
    DHS_similar_regions = DHS_hdf5['top_1000_similar_regions'][:]

    TF_hdf5 = h5py.File(TF_hdf5_path, 'r')
    TF_scaled_vals = TF_hdf5['top_1000_scaled_vals'][:]
    TF_similar_regions = TF_hdf5['top_1000_similar_regions'][:]
    
    epilogos_hdf5 = h5py.File(epilogos_hdf5_path, 'r')
    epilogos_scaled_vals = epilogos_hdf5['top_1000_scaled_vals'][:]
    epilogos_similar_regions = epilogos_hdf5['top_1000_similar_regions'][:]

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
    func = partial(count_terms, descriptor_terms, DHS_scaled_vals, DHS_similar_regions, TF_scaled_vals, TF_similar_regions, epilogos_scaled_vals, epilogos_similar_regions, mapped_genes)
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
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], int(sys.argv[8]))
