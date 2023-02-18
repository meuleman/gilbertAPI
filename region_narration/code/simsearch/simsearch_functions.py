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

    factor_order = np.append(DHS_factors, epi_factors)[target_factors]
    metric_order = np.append(np.array(['DHSs'] * len(DHS_factors)),
                             np.array(['chromatin states'] *
                                      len(epi_factors)))[target_factors]
    colors = np.append(DHS_color, epi_color)[target_factors]

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
