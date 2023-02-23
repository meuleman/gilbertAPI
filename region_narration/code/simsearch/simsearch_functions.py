import numpy as np
from datetime import date
import h5py
import os


def roiURL(roi, top_regions, link_type='dual'):
    today = date.today()
    roiurl_name = ('region_narration_simsearch_regions_roi_' + str(roi) +
                   '_scale_25kb_CT' + today.strftime("%Y%m%d") + '.bed')
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


def simsearch(DHS_per_f, epi_per_f, roi, target_factors):
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

    # convolution
    roi_percents = percents[roi, target_factors]
    roi_conv_all = percents[:, target_factors] * roi_percents

    # fusion
    roi_sim_all = 1 - np.prod(1 - roi_conv_all, axis=1)

    # most similar regions
    top_regions = np.argsort(-roi_sim_all)[:101]
    top_regions = np.append(roi, top_regions[top_regions != roi][:100])

    return (top_regions, roi_sim_all[top_regions], roi_conv_all[top_regions],
            percents[top_regions], factor_order, metric_order, colors)


def narrate_sim(top_regions, roi_per, factor_order, metric_order,
                target_factors, conversion_dict, metriccolors, thresh=0.9):

    query = top_regions[0]
    # convolution scores for all factors
    factor_scores = roi_per[0] * roi_per  # > sqrt of thresh for query region

    sentences = []
    for region, factor_score in zip(top_regions, factor_scores):
        # factors that rank high between both the query region and similar
        # region
        key_factors = np.where(factor_score > thresh)[0]
        if len(key_factors) > 0:
            if region == query:
                sentence = '<div>The query region ranks highly for '
            else:
                sentence = ('<div>This region is highly similar to the ' +
                            'query region in terms of ')
            # factors/metric to highlight
            key_factor_names = factor_order[key_factors]
            key_metric_names = metric_order[key_factors]
            m_sections = []
            # one metric at a time (unordered)
            for m in np.unique(key_metric_names):
                m_inds = np.where(key_metric_names == m)[0]
                # factors for current metric
                m_factors = key_factor_names[m_inds]
                # to bold target factors
                target_factor_inds = np.where(
                    np.isin(m_factors,
                            factor_order[target_factors])
                )[0]
                # convert to readable names
                m_factors = list(map(conversion_dict.get, m_factors))
                # bold target factors
                for t in target_factor_inds:
                    m_factors[t] = (metriccolors['BOLD'] + m_factors[t] +
                                    metriccolors['UNBOLD'])
                # metric color
                factors_joined = metriccolors[m]
                # join factors
                if len(m_factors) > 2:
                    factors_joined += ', '.join(m_factors[:-1])
                    factors_joined += ', and ' + m_factors[-1]
                elif len(m_factors) == 2:
                    factors_joined += ' and '.join(m_factors)
                else:
                    factors_joined += m_factors[0]

                factors_joined += ' ' + m + metriccolors['STANDARD']
                m_sections.append(factors_joined)
            # join metrics
            if len(m_sections) > 2:
                metrics_joined = ', '.join(m_sections[:-1])
                metrics_joined += ', and ' + m_sections[-1]
            if len(m_sections) == 2:
                metrics_joined = ' and '.join(m_sections)
            else:
                metrics_joined = m_sections[0]

            sentence += metrics_joined + '.</div>'

        else:
            sentence = '<div></div>'
        sentences.append(sentence)
    return sentences
