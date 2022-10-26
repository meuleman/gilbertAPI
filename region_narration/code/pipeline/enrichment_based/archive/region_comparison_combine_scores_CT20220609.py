import numpy as np
import sys
import os
import time
from datetime import date

def combine_scores(DHS_count_rank_vals, DHS_count_similar_regions, TF_count_rank_vals, TF_count_similar_regions, epilogos_score_rank_vals, epilogos_score_similar_regions):
    #max ind for combined score array
    max_ind = np.max([DHS_count_similar_regions, TF_count_similar_regions, epilogos_score_similar_regions])
    #DHS scores
    DHS_arr = np.zeros(max_ind + 1)
    DHS_arr[DHS_count_similar_regions] = DHS_count_rank_vals
    #epilogos scores
    epilogos_arr = np.zeros(max_ind + 1)
    epilogos_arr[epilogos_score_similar_regions] = epilogos_score_rank_vals
    #TF motif scores
    TF_arr = np.zeros(max_ind + 1)
    TF_arr[TF_count_similar_regions] = TF_count_rank_vals

    #scoring all metrics
    score_arr = 1 - ((1 - epilogos_arr) * (1 - DHS_arr) * (1 - TF_arr))

    return score_arr

def region_return(score_arr, roi, regions_path):
    #top 100 region inds and scores
    top_region_inds = np.append(roi, np.argsort(-score_arr)[:100])
    top_region_scores = np.append(1, -np.sort(-score_arr)[:100])
    #find genomic locations
    f = open(regions_path, 'r')
    regions = f.readlines()
    top_regions = []
    for ind in top_region_inds:
        top_regions.append(regions[ind].strip().split())
    f.close()
    top_regions = np.array(top_regions)
    top_regions = np.hstack((top_regions, top_region_scores.reshape((len(top_region_scores), 1))))

    link_path, bed_path_name = roi_url(roi, top_regions)
    print(link_path, flush=True)    

    return bed_path_name

def roi_url(roi, top_regions):
    today = date.today()
    roi_url_name = 'region_comparison_recommended_regions_roi_'+str(roi)+'_scale_250kb_CT' + today.strftime("%Y%m%d") + '.bed'
    bed_path_name = os.path.join('/home/ctrader/public_html/', roi_url_name)
    np.savetxt(bed_path_name, top_regions, delimiter='\t', fmt='%s')
    link_path = 'https://epilogos.altius.org:3001/?application=viewer&sampleSet=vA&mode=single&genome=hg19&model=15&complexity=KL&group=all&chrLeft=chr1&chrRight=chr1&start=8000&stop=35000&roiURL=https%3A%2F%2Fencode%3Acollabor8%40resources.altius.org%2F~ctrader%2F' + roi_url_name    
    return link_path, bed_path_name

def main(region_of_interest_path, DHS_npz_path, TF_npz_path, epilogos_npz_path, regions_path):
    t = time.time()
    #find the region of interest
    f = open(region_of_interest_path, 'r')
    roi = int(f.readline())

    ###Load rank and ind arrays
    #DHS counts
    DHS_count_rank_vals = np.load(DHS_npz_path)['top_1000_rank_vals'][roi]
    DHS_count_similar_regions = np.load(DHS_npz_path)['top_1000_similar_regions'][roi]

    #TF motif counts
    TF_count_rank_vals = np.load(TF_npz_path)['top_1000_rank_vals'][roi]
    TF_count_similar_regions = np.load(TF_npz_path)['top_1000_similar_regions'][roi]

    #epilogos score sum and max
    epilogos_score_rank_vals = np.load(epilogos_npz_path)['top_1000_rank_vals'][roi]
    epilogos_score_similar_regions = np.load(epilogos_npz_path)['top_1000_similar_regions'][roi]

    ###Combine scores 
    score_arr = combine_scores(DHS_count_rank_vals, DHS_count_similar_regions, TF_count_rank_vals, TF_count_similar_regions, epilogos_score_rank_vals, epilogos_score_similar_regions)

    ###Return regions
    bed_path_name = region_return(score_arr, roi, regions_path)
    print(bed_path_name)
    print(time.time()-t)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
