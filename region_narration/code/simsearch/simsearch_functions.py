import numpy as np
from datetime import date
import h5py
import os
import scipy.stats as stats


def roiURL(roi, top_regions, detail_level, target_factors_str,
           link_type='dual'):
    today = date.today()
    roiurl_name = ('region_narration_simsearch_regions_roi_' + str(roi) +
                   '_detail_level_' + str(detail_level) + '_factors_' +
                   target_factors_str + '_scale_25kb_CT' +
                   today.strftime("%Y%m%d") + '.bed')
    bed_path_name = os.path.join('/home/ctrader/public_html/', roiurl_name)
    np.savetxt(bed_path_name, top_regions, delimiter='\t', fmt='%s')
    link_path1 = ('https://epilogos.altius.org:1234/?application=viewer&' +
                  'sampleSet=vC&mode=single&genome=hg38&model=18&complex' +
                  'ity=KL&group=all&chrLeft=chr1&chrRight=chr1&start=800' +
                  '0&stop=35000&roiURL=https%3A%2F%2Fencode%3Acollabor8%' +
                  '40resources.altius.org%2F~ctrader%2F' + roiurl_name)
    link_path2 = ('https://index.altius.org:1234/?application=viewer&mod' +
                  'e=pd&chrLeft=chr11&chrRight=chr11&start=47329800&stop' +
                  '=47356800&roiURL=https%3A%2F%2Fencode%3Acollabor8%40r' +
                  'esources.altius.org%2F~ctrader%2F' + roiurl_name)
    link_path_dual = ('https://index.altius.org:3000/?application=viewer' +
                      '&mode=cf&&roiURL=https%3A%2F%2Fencode%3Acollabor8' +
                      '%40resources.altius.org%2F~ctrader%2F' + roiurl_name)
    if link_type == 'dual':
        return link_path_dual
    else:
        return link_path1, link_path2


def load_per(per_f):
    hdf5_f = h5py.File(per_f, 'r')
    percentile = hdf5_f['percentile'][:]
    factors = np.char.decode(hdf5_f['factors'][:])
    color = np.char.decode(hdf5_f['colors'][:])
    hdf5_f.close()
    return percentile, factors, color


def top_region_sig(percents, top_regions):
    # mean and std of factor percentiles in all regions
    all_mean = np.mean(percents, axis=0)
    all_std = np.std(percents, axis=0)

    # mean of factor percentiles in most similar regions
    top_region_mean = np.mean(percents[top_regions[:, 1:]], axis=1)

    # Z score
    Z = (top_region_mean - all_mean) / (all_std / np.sqrt(top_regions.shape[1] - 1))

    # p-value (one sided)
    pval = stats.norm.sf(Z)

    return pval, Z


def simsearch(percents, roi, included_factors, num_similar_regions):

        # factors the query region rank highest in
    query_factor_order = np.argsort(-percents[roi])

    # factors to search over
    target_factors = query_factor_order[
        np.isin(query_factor_order, included_factors)
    ]

    # convolution
    roi_percents = percents[roi, target_factors]
    roi_conv_all = percents[:, target_factors] * roi_percents

    # fusion
    roi_sim_all = []
    for dl in range(1, len(included_factors) + 1):
        roi_sim = 1 - np.prod(1 - roi_conv_all[:, :dl], axis=1)
        roi_sim_all.append(roi_sim)
    roi_sim_all = np.array(roi_sim_all)

    # most similar regions
    top_regions = np.zeros(
        (len(included_factors), num_similar_regions + 1)
        , dtype=int
    )
    top_regions_pre = np.argsort(
        -roi_sim_all, axis=1
    )[:, :(num_similar_regions + 1)]

    for i in range(len(included_factors)):
        top_regions[i] = np.append(
            roi,
            top_regions_pre[i][
                top_regions_pre[i] != roi
            ][:num_similar_regions]
        )

    # factor percentile significance
    pval, Z = top_region_sig(percents, top_regions)

    rows = (
        np.ones(top_regions.shape, dtype=int) *
        np.arange(top_regions.shape[0]).reshape(-1, 1)
    )

    return (top_regions, roi_sim_all[rows, top_regions],
            percents[top_regions], query_factor_order,
            target_factors, pval, Z)


def narrate_detail_level(top_regions, roi_per, factor_order, metric_order, included_factors, detail_level, conversion_dict, metriccolors, thresh=0.9):
    query = top_regions[0]
    hits = top_regions[1:]
    factor_scores = roi_per[0] * roi_per  # > sqrt of thresh for query region

    sentences = []
    for region, factor_score in zip(top_regions, factor_scores):
        factor_score_sort = np.argsort(-factor_score)
        key_factors = factor_score_sort[
            ((factor_score[factor_score_sort] > thresh) &
             (np.isin(factor_score_sort, included_factors)))
        ][:detail_level]

        if len(key_factors) > 0:
            if region == query:
                sentence = '<div>This region is characterized by '
            else:
                sentence = '<div>This region is highly similar to the query region in terms of '
            key_factor_names = factor_order[key_factors]
            key_factor_names = key_factor_names.astype(object)
            for i in range(len(key_factor_names)):
                 key_factor_names[i] = conversion_dict[key_factor_names[i]]
            key_metric_names = metric_order[key_factors]

            string_factors = []
            for factor, metric in zip(key_factor_names, key_metric_names):
                if metric == 'chromatin states':
                    string_factor = metriccolors['chromatin states']
                elif metric == 'DHSs':
                    string_factor = metriccolors['DHSs']
                string_factor += factor + ' ' + metric + metriccolors['STANDARD']
                string_factors.append(string_factor)

            if len(string_factors) > 2:
                sentence += ', '.join(string_factors[:-1])
                sentence += ', and ' + string_factors[-1]
            elif len(string_factors) == 2:
                sentence += ' and '.join(string_factors)
            else:
                sentence += string_factors[0]
            sentence += '.</div>'
        else:
            sentence = '<div></div>'
        sentences.append(sentence)
    return sentences
