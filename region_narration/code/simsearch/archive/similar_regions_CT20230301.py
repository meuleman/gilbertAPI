import numpy as np
import sys
import tabix
import json
from simsearch_functions import simsearch, narrate_detail_level


def main(DHS_per_f, epi_per_f, gz_file_path, regions_path,
         genomic_location, included_factors_str, detail_level):

    if str.isalpha(str(genomic_location)[0]):
        # load roi
        gz_tb = tabix.open(gz_file_path)
        roi = np.array(
            list(gz_tb.querys(genomic_location))
        )[:, 3].astype(int)[0]
    else:
        roi = int(genomic_location)

    # factors the user cares about
    included_factors = np.int_(
        included_factors_str.strip('[').strip(']').split(',')
    )

    # find most similar region inds
    (top_regions, roi_sim, roi_conv,
     roi_per, factor_order, metric_order,
     colors, query_factor_order,
     target_factors, pval,
     Z) = simsearch(DHS_per_f, epi_per_f,
                    roi, included_factors,
                    detail_level)

    # genomic coordinates for each region
    f = open(regions_path, 'r')
    regions = []
    for line in f.readlines():
        line = line.strip().split('\t')
        regions.append(line[0] + ':' + line[1] + '-' + line[2])
    f.close()

    top_region_coords = np.array(regions)[top_regions]

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

    # narrate regions
    sentences = narrate_detail_level(top_regions, roi_per,
                                     factor_order, metric_order,
                                     included_factors, detail_level, conversion_dict,
                                     metriccolors)

    # reorder based on query region
    # roi_conv = roi_conv[:, query_factor_order]
    roi_per = roi_per[:, query_factor_order]
    factor_order = factor_order[query_factor_order]
    metric_order = metric_order[query_factor_order]
    colors = colors[query_factor_order]
    pval = pval[query_factor_order]
    Z = Z[query_factor_order]

    # factor_order = factor_order.astype(object)
    for i in range(len(factor_order)):
        factor_order[i] = factor_order[i].replace('_', '/')
    #     factor_order[i] = metriccolors[metric_order[i]] + factor_order[i] + metriccolors['STANDARD']

    # for the roi, return results as json array
    # each element is an object of scores for each metric
    return_arr = []
    for rank in range(101):
        region_dict = {
            'rank': rank,
            'coordinates': top_region_coords[rank],
            'region_ind': int(top_regions[rank]),
            'score': float(roi_sim[rank]),
            # 'factor_conv': list(roi_conv[rank]),
            'percentiles': list(roi_per[rank]),
            'factors': list(factor_order),
            'metrics': list(metric_order),
            'factor_colors': list(colors),
            'sentence': sentences[rank],
        }
        if rank == 0:
            region_dict['query_factor_order'] = query_factor_order.tolist()
            region_dict['target_factors'] = target_factors.tolist()
            region_dict['top_region_factor_Z_score'] = list(Z)
            region_dict['top_region_factor_pval'] = list(pval)
        return_arr.append(region_dict)

    json_return = json.dumps(return_arr)
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
         sys.argv[6], int(sys.argv[7]))
