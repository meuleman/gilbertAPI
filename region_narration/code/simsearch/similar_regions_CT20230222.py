import numpy as np
import sys
import tabix
import json
from simsearch_functions import simsearch, narrate_sim


def main(DHS_per_f, epi_per_f, gz_file_path,
         regions_path, genomic_location, target_factors):

    # load roi
    gz_tb = tabix.open(gz_file_path)
    roi = np.array(list(gz_tb.querys(genomic_location)))[:, 3].astype(int)[0]

    # find most similar region inds
    (top_regions, roi_sim, roi_conv,
     roi_per, factor_order, metric_order,
     colors) = simsearch(
         DHS_per_f, epi_per_f, roi, target_factors
     )

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
    sentences = narrate_sim(top_regions, roi_per,
                            factor_order, metric_order,
                            target_factors, conversion_dict,
                            metriccolors)



    # for the roi, return results as json array
    # each element is an object of scores for each metric
    return_arr = []
    for rank in range(101):
        region_dict = {
            'rank': rank,
            'coordinates': top_region_coords[rank],
            'region_ind': int(top_regions[rank]),
            'score': float(roi_sim[rank]),
            'factor_conv': list(roi_conv[rank]),
            'percentiles': list(roi_per[rank]),
            'factors': list(factor_order),
            'metrics': list(metric_order),
            'factor_colors': list(colors),
            'sentence': sentences[rank]
        }
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
         np.int_(sys.argv[6].strip('[').strip(']').split(',')))
