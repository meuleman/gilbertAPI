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

def make_sentence(top_region_inds, DHS_arr_top_100, TF_arr_top_100, epilogos_arr_top_100, string_arr_top_100, DHS_roi_rank, TF_roi_rank,
                  epilogos_roi_mean_rank, epilogos_roi_max_rank, DHS_component_order, TF_component_order, epilogos_component_order, string_comp_order, 
                  DHSfact2descr, TFfact2descr, epilogosfact2descr, all_arr_ori, all_arr_dest, DHS_ooe, TF_ooe, epilogos_ooe, geneID2name):
    
    topfactor2region = {}
    #find top comps for each metric/recommended region
    #DHS scores
    for i in range(len(DHS_arr_top_100)):
        if DHS_arr_top_100[i] > 0:
            DHS_i_rank = rankdata(-DHS_ooe[top_region_inds[i+1]])
            top_comp = np.argsort((DHS_roi_rank+DHS_i_rank)/2)[0]
            f = 'DHS ' + DHSfact2descr[DHS_component_order[top_comp]]
            if f not in topfactor2region:
                topfactor2region[f] = np.array([i])
            else:
                topfactor2region[f] = np.append(topfactor2region[f], [i])
                
    #TF motif scores
    for i in range(len(TF_arr_top_100)):
        if TF_arr_top_100[i] > 0:
            TF_i_rank = rankdata(-TF_ooe[top_region_inds[i+1]])
            top_comp = np.argsort((TF_roi_rank+TF_i_rank)/2)[0]
            f = 'TF ' + TF_component_order[top_comp] + ': ' + TFfact2descr[TF_component_order[top_comp]]
            if f not in topfactor2region:
                topfactor2region[f] = np.array([i])
            else:
                topfactor2region[f] = np.append(topfactor2region[f], [i])

    #chromatin state scores
    for i in range(len(epilogos_arr_top_100)):
        if epilogos_arr_top_100[i] > 0:
            epilogos_i_mean_rank = rankdata(-epilogos_ooe[top_region_inds[i+1]][:18])
            epilogos_i_max_rank = rankdata(-epilogos_ooe[top_region_inds[i+1]][18:])
            top_mean_comp = np.argsort((epilogos_roi_mean_rank + epilogos_i_mean_rank)/2)[0]
            top_max_comp = np.argsort((epilogos_roi_max_rank + epilogos_i_max_rank)/2)[0]
            if top_mean_comp==top_max_comp:
                f = 'epilogos ' + epilogosfact2descr[epilogos_component_order[top_mean_comp]]
            else:
                f = 'epilogos ' + (epilogosfact2descr[epilogos_component_order[top_mean_comp]] + ', with ' + 
                                   epilogosfact2descr[epilogos_component_order[top_max_comp]] + ' peaks')
            if f not in topfactor2region:
                topfactor2region[f] = np.array([i])
            else:
                topfactor2region[f] = np.append(topfactor2region[f], [i])

    #STRING Scores
    for i in range(len(string_arr_top_100[0, :])):
        if np.sum(string_arr_top_100[:, i]) > 0:
            pos_string_inds = np.where(string_arr_top_100[:, i]>0)[0]
            for comp in pos_string_inds:
                ori_gene = geneID2name[all_arr_ori[comp, i]]
                f = 'Gene ' + string_comp_order[comp] + ' ' + ori_gene
                if f not in topfactor2region:
                    topfactor2region[f] = np.array([i])
                else:
                    topfactor2region[f] = np.append(topfactor2region[f], [i])
    
    #ordering facors by region counts
    factor_list = np.array([])
    region_counts_per_factor = np.array([], dtype=int)
    for factor in topfactor2region.keys():
        factor_list = np.append(factor_list, factor)
        region_counts_per_factor = np.append(region_counts_per_factor, len(topfactor2region[factor]))

    sort_order = np.argsort(-region_counts_per_factor)
    factor_list = factor_list[sort_order]
    region_counts_per_factor = region_counts_per_factor[sort_order]

    #finding factors representing top x percent of regions
    limiting_per = 95

    factor_counts_per_region = np.zeros(100, dtype=int)
    factor_list_top = np.array([])
    region_counts_per_factor_top = np.array([], dtype=int)
    for factor, num_regions in zip(factor_list, region_counts_per_factor):
        factor_counts_per_region[topfactor2region[factor]] += 1
        factor_list_top = np.append(factor_list_top, factor)
        region_counts_per_factor_top = np.append(region_counts_per_factor_top, num_regions)
        if np.count_nonzero(factor_counts_per_region) >= limiting_per:
            break
    
    ####
    #write a sentence based on 95%
    DHS_factors = []
    DHS_factor_counts = []
    TF_factors = []
    TF_factor_counts = []
    epi_factors = []
    epi_factor_counts = []
    string_factors = []
    string_factor_counts = []
    for factor, count in zip(factor_list_top, region_counts_per_factor_top):
        if factor.split()[0] == 'DHS':
            DHS_factors.append(' '.join(factor.split()[5:-1]))
            DHS_factor_counts.append(count)
        elif factor.split()[0] == 'TF':
            TF_factors.append(factor.split()[1][:-1])
            TF_factor_counts.append(count)
        elif factor.split()[0] == 'epilogos':
            epi_factors.append(' '.join(factor.split()[1:]))
            epi_factor_counts.append(count)
        elif factor.split()[0] == 'Gene':
            string_factors.append([factor.split()[-1], ' '.join(factor.split()[1:-1])])
            string_factor_counts.append(count)
    string_factors = np.array(string_factors)

    sentence = ''

    #Gene/STRING factors
    if len(string_factors)>0:
        for gene in np.unique(string_factors[:, 0]):
            gene_inds = np.where(string_factors[:, 0]==gene)[0]
            sentence += 'Gene ' + gene + ' (with '
            for ind in gene_inds[:-1]:
                sentence += string_factors[ind, 1] + ' in ' + str(string_factor_counts[ind]) + '% of regions, '
            sentence += string_factors[gene_inds[-1], 1] + ' in ' + str(string_factor_counts[gene_inds[-1]]) + '% of regions)'

        if ((len(DHS_factors)>0) | (len(TF_factors)>0) | (len(epi_factors)>0)):
            sentence += ', '
        else:
            sentence += '.'

    #DHS factors
    if len(DHS_factors)>0:
        sentence += 'DNase I accessibility (in '
        for factor, count in zip(DHS_factors[:-1], DHS_factor_counts[:-1]):
            sentence += factor + ' samples in ' + str(count) + '% of regions, '
        sentence += DHS_factors[-1] + ' samples in ' + str(DHS_factor_counts[-1]) + '% of regions)'
        if ((len(TF_factors)>0) | (len(epi_factors)>0)):
            sentence += ', '
        else:
            sentence += '.'

    #TF factors
    if len(TF_factors)>0:
        sentence += 'Transcription Factor motifs ('
        for factor, count in zip(TF_factors[:-1], TF_factor_counts[:-1]):
            sentence += factor + ' in ' + str(count) + '% of regions, '
        sentence += TF_factors[-1] + ' in ' + str(TF_factor_counts[-1]) + '% of regions)'
        if ((len(epi_factors)>0)):
            sentence += ', '
        else:
            sentence += '.'

    #Chromatin States
    if len(epi_factors)>0:
        sentence += 'Chromatin States ('
        for factor, count in zip(epi_factors[:-1], epi_factor_counts[:-1]):
            sentence += factor + ' in ' + str(count) + '% of regions, '
        sentence += epi_factors[-1] + ' in ' + str(epi_factor_counts[-1]) + '% of regions)'
        sentence += '.'

    return sentence

def main(gz_file_path, genomic_location, DHS_hdf5_path, TF_hdf5_path, epilogos_hdf5_path, coex_hdf5_path, cooccur_hdf5_path, database_hdf5_path, exper_hdf5_path,
         fusion_hdf5_path, homology_hdf5_path, neigh_hdf5_path, textmin_hdf5_path, regions_path, DHS_ooe_path, TF_ooe_path, epilogos_ooe_path,
         DHS_mapped_path, TF_mapped_path, epilogos_mapped_path, DHSfact2descr_path,  TFfact2descr_path, epifact2descr_path, geneID2name_path, out_bed_path, sentence_out):

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
    string_comp_order = ['Co-expression', 'Co-occurence', 'Co-annotation', 
                         'Interactions in vitro', 'Fusion in Other Species', 'Common Ancestory', 
                         'being Nearby', 'Co-citation']

    #roi comps
    DHS_roi_rank = rankdata(-DHS_ooe[top_region_inds[:, 0]], axis=1)
    TF_roi_rank = rankdata(-TF_ooe[top_region_inds[:, 0]], axis=1)
    epilogos_roi_mean_rank = rankdata(-epilogos_ooe[top_region_inds[:, 0]][:, :18], axis=1)
    epilogos_roi_max_rank = rankdata(-epilogos_ooe[top_region_inds[:, 0]][:, 18:], axis=1)

    #find genomic locations
    f = open(regions_path, 'r')
    regions = f.readlines()
    regions_for_sentences = []
    for ind in roi:
        regions_for_sentences.append(regions[ind].strip().split())
    f.close()

    #generate sentences
    file = open(sentence_out, 'w')
    for r in range(len(roi)):
        sentence = make_sentence(top_region_inds[r], DHS_arr_top_100[r], TF_arr_top_100[r], epilogos_arr_top_100[r], string_arr_top_100[:, r], DHS_roi_rank[r], TF_roi_rank[r],
                                 epilogos_roi_mean_rank[r], epilogos_roi_max_rank[r], DHS_component_order, TF_component_order, epilogos_component_order, string_comp_order, 
                                 DHSfact2descr, TFfact2descr, epilogosfact2descr, all_arr_ori[:, r], all_arr_dest[:, r], DHS_ooe, TF_ooe, epilogos_ooe, geneID2name)
        print(roi[r], r)
        print(sentence, flush=True)
        file.write(regions_for_sentences[r][0] + '\t' + regions_for_sentences[r][1] + '\t' + 
                   regions_for_sentences[r][2] + '\t' + str(roi[r]) + '\t' + sentence + '\n')

    file.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12],
         sys.argv[13], sys.argv[14], sys.argv[15], sys.argv[16], sys.argv[17], sys.argv[18], sys.argv[19], sys.argv[20], sys.argv[21], sys.argv[22], sys.argv[23], 
         sys.argv[24], sys.argv[25], sys.argv[26])
