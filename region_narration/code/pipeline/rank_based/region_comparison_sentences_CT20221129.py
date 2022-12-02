import h5py
import numpy as np
"""
Functions for sentence generation.
"""



def load_rank_path(hdf5_path, roi):
    hdf5_file = h5py.File(hdf5_path, 'r')
    rank_arr = hdf5_file['ranked_regions'][roi]
    max_ranks = hdf5_file['max_ranks'][:]
    hdf5_file.close()
    return rank_arr, max_ranks

def get_factor_order(path):
    if path[-5:]=='.hdf5':
        hdf5_file = h5py.File(path, 'r')
        factor_order = np.char.decode(hdf5_file['comp_order'][:])
        hdf5_file.close()
    else:
        f = open(path, 'r')
        factor_order = np.array(f.readline().strip().split('\t'))[3:]
        f.close()
    return factor_order

def regional_sentence(rank_arrs, factor_arrs, max_rank_arrs, metric_name_arrs,
                      conversion_dicts, max_rank=100, multiple_rank_thresh=5,
                      max_factors_per_metric=3):

    num_metrics = len(rank_arrs)
    num_regions = len(rank_arrs[0])

    #establish sentences
    sentences = np.array([''] * num_regions, dtype=object)
    sentence_start = 'This region ranks in the top '
    sentences[:] = sentence_start

    sentence_section_arr = []
    sentence_section_best_ranks = []

    backup_rank = -np.ones(num_regions)
    backup_factor = np.array([''] * num_regions, dtype=object)
    backup_metric = np.array([''] * num_regions, dtype=object)

    #ranks for each metric
    for (rank_arr, factor_arr,
         max_ranks, metric_name,
         conversion_dict, metric_i) in zip(rank_arrs, factor_arrs,
                                           max_rank_arrs,
                                           metric_name_arrs,
                                           conversion_dicts,
                                           range(num_metrics)[:-1]):

        #find percent of rank compared to max rank of each factor
        adjusted_ranks = rank_arr / (max_ranks + 1) #pseudocount
        #find the top ranking factor
        rank_sort_order = np.argsort(-adjusted_ranks, axis=1)
        ranks_sorted = -np.sort(-adjusted_ranks, axis=1)

        #int of rank percentage
        #ranks_percent = np.ceil((1 - ranks_sorted) * 100).astype(int)
        ranks_percent = (1 - ranks_sorted) * 100
        ranks_percent = np.where(ranks_percent < 1, np.ceil(ranks_percent * 100) / 100, np.ceil(ranks_percent))
        #set top 0% ranks to top 1% ranks
    #     ranks_percent = np.where(ranks_percent==0, 1, ranks_percent) #not needed with pseudocount

        ###the ranks to include
        #include all values that rank in the top multiple_rank_thresh percent
        ranks_mask = ranks_percent <= multiple_rank_thresh
        #include top factor if rank is in the top max_rank percent
        ranks_mask[:, 0] = ((ranks_mask[:, 0]) | (ranks_percent[:, 0] <= max_rank))
        #number of factors to include for each region
        num_factors_per_region = np.count_nonzero(ranks_mask, axis=1)

        #find name of top factor
        factors_ranked = factor_arr[rank_sort_order]
        #convert to a readable version of name
        factors_ranked = np.array(list(map(conversion_dict.get, factors_ranked.flatten()))).reshape(factors_ranked.shape)

        #metric in sentence form
        metric_arr = np.char.add(' ', np.array([metric_name] * num_regions))

        #backup for regions
        better_backup_mask = ranks_sorted[:, 0] > backup_rank
        backup_rank[better_backup_mask] = ranks_sorted[:, 0][better_backup_mask]
        backup_factor[better_backup_mask] = factors_ranked[:, 0][better_backup_mask]
        backup_metric[better_backup_mask] = metric_arr[better_backup_mask]

        #initializing the section of the sentence for current metric
        sentence_section = np.array([''] * num_regions)

        metric_ranks_per_region = []
        #for each factor in order of how they rank for each region, add to sentence section
        for factor in range(len(factor_arr))[:max_factors_per_metric]:
            if ((len(ranks_mask[:, factor][ranks_mask[:, factor]]) == 0) & (factor!=0)):
                break
            #mask for factor
            factor_mask = ~ranks_mask[:, factor]
            #percentage
            factor_rank_percent = ranks_percent[:, factor].astype(str)
            factor_rank_percent[factor_mask] = ''

            #saving the rank values for each region
            metric_ranks_per_region.append(np.where(factor_mask, np.nan, ranks_sorted[:, factor]))

            #filler between percent and factor name
            per_of = np.array(['% of '] * num_regions)
            per_of[factor_mask] = ''
            #factor name
            factor_name = factors_ranked[:, factor]
            factor_name[factor_mask] = ''

            #ending
            factor_ending = np.array([''] * num_regions, dtype=object)
            #for regions with more factors than max_factors_per_metric with ranks better than multiple_rank_thresh
            if (factor + 2) == max_factors_per_metric: #if there is 1 factor left
                if factor==0:
                    factor_ending[num_factors_per_region >= (factor + 2)] = ' and ' #>
                else:
                    factor_ending[num_factors_per_region >= (factor + 2)] = ', and ' #>
            elif (factor + 1) == max_factors_per_metric: #if there are no factors left
                factor_ending[num_factors_per_region >= (factor + 2)] = '' #>
            else: #many factors left
                factor_ending[num_factors_per_region > (factor + 2)] = ', '
                if factor==0:
                    factor_ending[num_factors_per_region==(factor + 2)] = ' and '
                else:
                    factor_ending[num_factors_per_region==(factor + 2)] = ', and '

    #         #for regions with 1 more factor left under the multiple_rank_thresh
    #         if (factor + 1) == max_factors_per_metric: #if there are no factors left
    #             factor_ending[num_factors_per_region==(factor + 2)] = ''
    #         else:
    #             factor_ending[num_factors_per_region==(factor + 2)] = ' and '

            #append factor info to sentence section
            sentence_section = np.char.add(sentence_section, factor_rank_percent)
            sentence_section = np.char.add(sentence_section, per_of)
            sentence_section = np.char.add(sentence_section, factor_name)
            sentence_section = np.char.add(sentence_section, factor_ending)

        #add the metric to the sentence section
        metric_arr[~ranks_mask[:, 0]] = ''
        sentence_section = np.char.add(sentence_section, metric_arr)

        sentence_section_arr.append(sentence_section)

        #best ranks for each region
        metric_ranks_per_region = np.array(metric_ranks_per_region).T
        sentence_section_best_ranks.append(np.nanmax(metric_ranks_per_region, axis=1))

    sentence_section_arr = np.array(sentence_section_arr).T
    sentence_section_best_ranks = np.array(sentence_section_best_ranks).T
    sorted_sections = np.argsort(-sentence_section_best_ranks, axis=1)

    #to know where to place '.', ',' and 'and'
    num_sections_with_ranks = np.count_nonzero(~np.isnan(sentence_section_best_ranks), axis=1)

    #add sections to sentences
    for section in range(num_metrics)[:-1]:
        sentences = np.char.add(sentences.astype(str),
                                sentence_section_arr[np.arange(len(sorted_sections)),
                                                     sorted_sections[:, section]])
        section_ending = np.array([', '] * num_regions, dtype=object)
        section_ending[num_sections_with_ranks==(section + 1)] = '.'
        section_ending[num_sections_with_ranks<(section + 1)] = ''
        section_ending[(num_sections_with_ranks==(section + 2)) &
                       (num_sections_with_ranks!=2)] = ', and '
        section_ending[(num_sections_with_ranks==(section + 2)) &
                       (num_sections_with_ranks==2)] = ' and '
        sentences = np.char.add(sentences.astype(str),
                                section_ending)

    #fill with backup ranks for regions that don't match filters for any metric
    need_backup_mask = sentences==sentence_start
    #rank
    backup_percent = np.ceil((1 - backup_rank) * 100).astype(int)
    sentences[need_backup_mask] = np.char.add(sentences[need_backup_mask].astype(str),
                                              backup_percent[need_backup_mask].astype(int).astype(str))
    #fill between rank and factor
    sentences[need_backup_mask] = np.char.add(sentences[need_backup_mask].astype(str),
                                              np.array(['% of '] * num_regions)[need_backup_mask])
    #factor
    sentences[need_backup_mask] = np.char.add(sentences[need_backup_mask].astype(str),
                                              backup_factor[need_backup_mask])
    #metric
    sentences[need_backup_mask] = np.char.add(sentences[need_backup_mask].astype(str),
                                              backup_metric[need_backup_mask])
    #.
    sentences[need_backup_mask] = np.char.add(sentences[need_backup_mask].astype(str),
                                              np.array(['.'] * num_regions)[need_backup_mask])
    return sentences

def make_sentences(DHS_hdf5_path, DHS_mapped_path, TF_hdf5_path, TF_mapped_path, epilogos_hdf5_path, epilogos_mapped_path, STRING_ranks_hdf5_path, roi):
    #load rank arrays
    DHS_rank_arr, DHS_max_ranks = load_rank_path(DHS_hdf5_path, roi)
    TF_rank_arr, TF_max_ranks = load_rank_path(TF_hdf5_path, roi)
    epi_rank_arr, epi_max_ranks = load_rank_path(epilogos_hdf5_path, roi)
    STRING_rank_arr, STRING_max_ranks = load_rank_path(STRING_ranks_hdf5_path, roi)
    STRING_rank_arr = STRING_rank_arr[:, -1].reshape(-1, 1)
    STRING_max_ranks = STRING_max_ranks[-1]

    #load factor orders
    DHS_factor_order = get_factor_order(DHS_mapped_path)
    TF_factor_order = get_factor_order(TF_mapped_path)
    epi_factor_order = get_factor_order(epilogos_mapped_path)
    STRING_factor_order = get_factor_order(STRING_ranks_hdf5_path)[[-1]]

    #conversion dictionaries
    epi_conversion_dict = {
        'EnhA1': 'active enhancer 1',
        'EnhA2': 'active enhancer 2',
        'EnhBiv': 'bivalent enhancer',
        'EnhG1': 'genic enhancer 1',
        'EnhG2': 'genic enhancer 2',
        'EnhWk': 'weak enhancer',
        'Het': 'heterochromatin',
        'Quies': 'quiescent',
        'ReprPC': 'repressed polycomb',
        'ReprPCWk': 'weak repressed polycomb',
        'TssA': 'active TSS',
        'TssBiv': 'bivalent TSS',
        'TssFlnk': 'flanking TSS',
        'TssFlnkD': 'downstream flanking TSS',
        'TssFlnkU': 'upstream flanking TSS',
        'Tx': 'strong transcription',
        'TxWk': 'weak transcription',
        'ZNF_Rpts': 'ZNF genes + repeats'
    }

    DHS_conversion_dict = {
        'Cancer_epithelial': 'cancer/epithelial',
        'Cardiac': 'cardiac',
        'Digestive': 'digestive',
        'Lymphoid': 'lymphoid',
        'Musculoskeletal': 'musculoskeletal',
        'Myeloid_erythroid': 'myeloid/erythroid',
        'Neural': 'neural',
        'Organdevel_renal': 'organ development/renal',
        'Placental_trophoblast': 'placental/trophoblast',
        'Primitive_embryonic': 'primitive/embryonic',
        'Pulmonarydevel': 'pulmonary development',
        'Renal_cancer': 'renal/cancer',
        'StromalA': 'stromal A',
        'StromalB': 'stromal B',
        'Tissueinvariant': 'tissue invariant',
        'Vascular_endothelial': 'vascular/endothelial'
    }

    TF_conversion_dict = {}
    for factor in TF_factor_order:
        TF_conversion_dict[factor] = factor

    STRING_conversion_dict = {
        'coexpression' : 'co-expression',
        'cooccurence' : 'co-evolution',
        'database' : 'co-annotation',
        'experiments' : 'in vitro',
        'fusion' : 'non-human fused',
        'homology' : 'common ancestry',
        'neighborhood' : 'neighboring',
        'textmining' : 'literature based',
        'combined_score': 'combined'
    }

    rank_arrs = [DHS_rank_arr, TF_rank_arr, epi_rank_arr[:, :18], epi_rank_arr[:, 18:], STRING_rank_arr]
    factor_arrs = [DHS_factor_order, TF_factor_order, epi_factor_order, epi_factor_order, STRING_factor_order]
    max_rank_arrs = [DHS_max_ranks, TF_max_ranks, epi_max_ranks[:18], epi_max_ranks[18:], STRING_max_ranks]
    metric_name_arrs = ['DHS components', 'TF motifs', 'chromatin state themes', 'chromatin state peaks', 'inter-gene associations']
    conversion_dicts = [DHS_conversion_dict, TF_conversion_dict, epi_conversion_dict, epi_conversion_dict, STRING_conversion_dict]

    sentences = regional_sentence(rank_arrs, factor_arrs, max_rank_arrs, metric_name_arrs,
                                  conversion_dicts, max_rank=10, multiple_rank_thresh=3,
                                  max_factors_per_metric=3)

    return sentences

#if __name__ == "__main__":
#    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13], sys.argv[14], sys.argv[15])
