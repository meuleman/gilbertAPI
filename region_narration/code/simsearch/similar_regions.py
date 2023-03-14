import numpy as np
import sys
import tabix
import json
from simsearch_functions import simsearch, narrate_detail_level, roiURL, load_per


def main(DHS_per_f, epi_per_f, gz_file_path, regions_path,
         genomic_location, included_factors_str, get_roiURL):

    if str.isalpha(str(genomic_location)[0]):
        # load roi
        gz_tb = tabix.open(gz_file_path)
        roi = np.array(
            list(gz_tb.querys(genomic_location))
        )[:, 3].astype(int)[0]
    else:
        roi = int(genomic_location)

    # load region percentiles
    DHS_per, DHS_factors, DHS_color = load_per(DHS_per_f)
    epi_per, epi_factors, epi_color = load_per(epi_per_f)

    # join data from each metric
    percents = np.hstack((DHS_per, epi_per))
    factor_order = np.append(DHS_factors, epi_factors)
    metric_order = np.append(np.array(['DHSs'] * len(DHS_factors)),
                             np.array(['chromatin states'] *
                                      len(epi_factors)))
    colors = np.append(DHS_color, epi_color)
    # conversion dictionary
    conversion_dict = {
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
        'ZNF_Rpts': 'ZNF genes + repeats',
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

    metriccolors = {
        'DHSs': '<span style="color:#0000CD">',
        'chromatin states': '<span style="color:#FF0000">',
        'STANDARD': '</span>',
        'BOLD': '<b>',
        'UNBOLD': '</b>'
    }

    # factors the user cares about
    included_factors = np.int_(
        included_factors_str.strip('[').strip(']').split(',')
    )

    # genomic coordinates for each region
    f = open(regions_path, 'r')
    regions = []
    for line in f.readlines():
        line = line.strip().split('\t')
        regions.append(line[0] + ':' + line[1] + '-' + line[2])
    f.close()

    # find most similar region inds
    (top_regions, roi_sim,
     roi_per, query_factor_order,
     target_factors, pval,
     Z) = simsearch(percents, roi, included_factors)

    # get coordinates
    top_region_coords = np.array(regions)[top_regions]

    # dual track roiURL
    roiURLs = []
    for dl in range(1, 35):
        if get_roiURL == 'True':
            roiURLs.append('roiURL not available')
            # roiURLs.append(roiURL(roi, top_region_coords[dl-1], dl,
            #                       '_'.join(target_factors.astype(str))))
        else:
            roiURLs.append('roiURL not available')

    # narrate regions
    sentences = []
    for i in range(34):
        detail_level = i + 1
        sentences_pre = narrate_detail_level(top_regions[i], roi_per[i],
                                             factor_order, metric_order,
                                             included_factors, detail_level, conversion_dict,
                                             metriccolors)
        sentences.append(sentences_pre)
    sentences = np.array(sentences)

    # reorder based on query region
    roi_per = roi_per[:, :, query_factor_order]
    return_factor_order = factor_order[query_factor_order]
    return_metric_order = metric_order[query_factor_order]
    return_colors = colors[query_factor_order]
    pval = pval[:, query_factor_order]
    Z = Z[:, query_factor_order]

    # return_factor_order = return_factor_order.astype(object)
    for i in range(len(return_factor_order)):
        return_factor_order[i] = return_factor_order[i].replace('_', '/')
    #     return_factor_order[i] = metriccolors[return_metric_order[i]] + return_factor_order[i] + metriccolors['STANDARD']

    return_arr_full = []
    for dl in range(34):
        # for the roi, return results as json array
        # each element is an object of scores for each metric
        return_arr = []
        for rank in range(101):
            region_dict = {
                'rank': rank,
                'coordinates': top_region_coords[dl, rank],
                'region_ind': int(top_regions[dl, rank]),
                'score': float(roi_sim[dl, rank]),
                'percentiles': list(roi_per[dl, rank]),
                'factors': list(return_factor_order),
                'metrics': list(return_metric_order),
                'factor_colors': list(return_colors),
                'sentence': sentences[dl, rank],
                'detail_level': dl + 1
            }
            if rank == 0:
                region_dict['query_factor_order'] = query_factor_order.tolist()
                region_dict['target_factors'] = target_factors.tolist()
                region_dict['top_region_factor_Z_score'] = list(Z[dl])
                region_dict['top_region_factor_pval'] = list(pval[dl])
                region_dict['roiURL'] = roiURLs[dl]
            return_arr.append(region_dict)

        return_arr_full.append(return_arr)

    json_return = json.dumps(return_arr_full)
    sys.stdout.write(json_return)

    # save output for dev mode
    # f = open('/Users/ctrader/Documents/AG/NarrativeGenomics/' +
    #          'region_narration/code/pipeline/region_similarity_react/' +
    #          'region-narration-app/src/' +
    #          'devmode_top_similar_regions_results_CT20221129.json', 'w')
    # f.write(json_return)
    # f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
         sys.argv[6], sys.argv[7])
