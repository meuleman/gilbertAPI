import numpy as np
import sys
import os
import time
from datetime import date
import h5py
import tabix

def combine_scores(DHS_count_scaled_vals, DHS_count_similar_regions, TF_count_scaled_vals, TF_count_similar_regions, epilogos_score_scaled_vals, epilogos_score_similar_regions):
    #max ind for combined score array
    max_ind = np.max([DHS_count_similar_regions, TF_count_similar_regions, epilogos_score_similar_regions])

    #DHS scores
    DHS_arr = np.zeros(max_ind + 1)
    DHS_arr[DHS_count_similar_regions] = DHS_count_scaled_vals

    #epilogos scores
    epilogos_arr = np.zeros(max_ind + 1)
    epilogos_arr[epilogos_score_similar_regions] = epilogos_score_scaled_vals

    #TF motif scores
    TF_arr = np.zeros(max_ind + 1)
    TF_arr[TF_count_similar_regions] = TF_count_scaled_vals

    #scoring all metrics
    score_arr = 1 - ((1 - epilogos_arr) * (1 - DHS_arr) * (1 - TF_arr))

    #labeling most contributing metrics for each region
    contributing_cosine_arr = np.argmax(np.vstack((DHS_arr, TF_arr, epilogos_arr)), axis=0)

    return score_arr, contributing_cosine_arr

def region_return(score_arr, contributing_cosine_arr, roi, regions_path):
    #top 100 region inds and scores
    rec_inds = np.argsort(-score_arr)[:100]
    top_region_inds = np.append(roi, rec_inds)
    top_region_scores = np.append(1, score_arr[rec_inds])

    #most contributing metrics to the top 100 scoring regions
    top_region_contributions = np.append('ROI', contributing_cosine_arr[rec_inds])
    top_region_contributions = np.where(top_region_contributions=='0', 'DHS', top_region_contributions)
    top_region_contributions = np.where(top_region_contributions=='1', 'TF', top_region_contributions)
    top_region_contributions = np.where(top_region_contributions=='2', 'epilogos', top_region_contributions)
    
    #combining scores and contributions
    score_contribution_comb = np.char.add(np.char.add(np.round(top_region_scores, 4).astype(str), 
                                                      np.array([' - ']*len(top_region_scores))), 
                                          top_region_contributions)

    #find genomic locations
    f = open(regions_path, 'r')
    regions = f.readlines()
    top_regions = []
    for ind in top_region_inds:
        top_regions.append(regions[ind].strip().split())
    f.close()
    top_regions = np.array(top_regions)
    top_regions = np.hstack((top_regions, score_contribution_comb.reshape((len(score_contribution_comb), 1))))

    link_path, bed_path_name = roi_url(roi, top_regions)
    print(link_path, flush=True)    

    return bed_path_name

def roi_url(roi, top_regions):
    today = date.today()
    roi_url_name = 'region_comparison_recommended_regions_roi_'+str(roi)+'_scale_250kb_CT' + today.strftime("%Y%m%d") + '.bed'
    bed_path_name = os.path.join('/home/ctrader/public_html/', roi_url_name)
    np.savetxt(bed_path_name, top_regions, delimiter='\t', fmt='%s')
    link_path = 'https://epilogos.altius.org:3001/?application=viewer&sampleSet=vA&mode=single&genome=hg38&model=18&complexity=KL&group=all&chrLeft=chr1&chrRight=chr1&start=8000&stop=35000&roiURL=https%3A%2F%2Fencode%3Acollabor8%40resources.altius.org%2F~ctrader%2F' + roi_url_name    
    return link_path, bed_path_name

def main(gz_file_path, genomic_location, DHS_hdf5_path, TF_hdf5_path, epilogos_hdf5_path, regions_path):
    t = time.time()

    #find roi
    gz_tb = tabix.open(gz_file_path)
    roi = int(list(gz_tb.querys(genomic_location))[0][3])

    print('roi = ' + str(roi), flush=True)

    ###Load rank and ind arrays
    #DHS counts
    DHS_hdf5 = h5py.File(DHS_hdf5_path, 'r')
    DHS_count_scaled_vals = DHS_hdf5['top_1000_scaled_vals'][roi]
    DHS_count_similar_regions = DHS_hdf5['top_1000_similar_regions'][roi]
    DHS_hdf5.close()

    #TF motif counts
    TF_hdf5 = h5py.File(TF_hdf5_path, 'r')
    TF_count_scaled_vals = TF_hdf5['top_1000_scaled_vals'][roi]
    TF_count_similar_regions = TF_hdf5['top_1000_similar_regions'][roi]
    TF_hdf5.close()

    #epilogos score sum and max
    epilogos_hdf5 = h5py.File(epilogos_hdf5_path, 'r')
    epilogos_score_scaled_vals = epilogos_hdf5['top_1000_scaled_vals'][roi]
    epilogos_score_similar_regions = epilogos_hdf5['top_1000_similar_regions'][roi]
    epilogos_hdf5.close()

    ###Combine scores 
    score_arr, contributing_cosine_arr = combine_scores(DHS_count_scaled_vals, DHS_count_similar_regions, TF_count_scaled_vals, TF_count_similar_regions, epilogos_score_scaled_vals, epilogos_score_similar_regions)

    ###Return regions
    bed_path_name = region_return(score_arr, contributing_cosine_arr, roi, regions_path)
    print(bed_path_name)
    print(time.time()-t)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
