import numpy as np
import sys
import os
import time
from datetime import date
import h5py
import tabix
from scipy.stats import rankdata

def combine_scores(DHS_count_scaled_vals, DHS_count_similar_regions, TF_count_scaled_vals, TF_count_similar_regions, epilogos_score_scaled_vals,
                   epilogos_score_similar_regions, coex_scaled_vals, coex_similar_regions, cooccur_scaled_vals, cooccur_similar_regions, 
		   database_scaled_vals, database_similar_regions, exper_scaled_vals, exper_similar_regions,fusion_scaled_vals, fusion_similar_regions, 
		   homology_scaled_vals, homology_similar_regions, neigh_scaled_vals, neigh_similar_regions, textmin_scaled_vals, textmin_similar_regions):
    #max ind for combined score array
    max_ind = np.max([DHS_count_similar_regions, TF_count_similar_regions, epilogos_score_similar_regions, coex_similar_regions,
		     cooccur_similar_regions, database_similar_regions, exper_similar_regions, fusion_similar_regions, homology_similar_regions,
		     neigh_similar_regions, textmin_similar_regions])

    #DHS scores
    DHS_arr = np.zeros(max_ind + 1)
    DHS_arr[DHS_count_similar_regions] = DHS_count_scaled_vals
    #epilogos scores
    epilogos_arr = np.zeros(max_ind + 1)
    epilogos_arr[epilogos_score_similar_regions] = epilogos_score_scaled_vals
    #TF motif scores
    TF_arr = np.zeros(max_ind + 1)
    TF_arr[TF_count_similar_regions] = TF_count_scaled_vals
    #coexpression scores
    coex_arr = np.zeros(max_ind + 1)
    coex_arr[coex_similar_regions] = coex_scaled_vals
    #cooccurrence scores
    cooccur_arr = np.zeros(max_ind + 1)
    cooccur_arr[cooccur_similar_regions] = cooccur_scaled_vals
    #database scores
    database_arr = np.zeros(max_ind + 1)
    database_arr[database_similar_regions] = database_scaled_vals
    #experiment scores
    exper_arr = np.zeros(max_ind + 1)
    exper_arr[exper_similar_regions] = exper_scaled_vals
    #fusion scores
    fusion_arr = np.zeros(max_ind + 1)
    fusion_arr[fusion_similar_regions] = fusion_scaled_vals
    #homology scores
    homology_arr = np.zeros(max_ind + 1)
    homology_arr[homology_similar_regions] = homology_scaled_vals
    #neighborhood scores
    neigh_arr = np.zeros(max_ind + 1)
    neigh_arr[neigh_similar_regions] = neigh_scaled_vals
    #text mining scores
    textmin_arr = np.zeros(max_ind + 1)
    textmin_arr[textmin_similar_regions] = textmin_scaled_vals

    #scoring all metrics
    score_arr = 1 - ((1 - epilogos_arr) * (1 - DHS_arr) * (1 - TF_arr) * (1 - coex_arr) * (1 - cooccur_arr) * (1 - database_arr) * (1 - exper_arr) * 
    		     (1 - fusion_arr) * (1 - homology_arr) * (1 - neigh_arr) * (1 - textmin_arr))

    #labeling most contributing metrics for each region
    contributing_cosine_arr = np.argmax(np.vstack((DHS_arr, TF_arr, epilogos_arr, coex_arr, cooccur_arr, database_arr, exper_arr, fusion_arr, 
						   homology_arr, neigh_arr, textmin_arr)), axis=0)

    return score_arr, contributing_cosine_arr

def region_return(score_arr, contributing_cosine_arr, roi, regions_path, DHS_ooe_path, TF_ooe_path, epilogos_ooe_path, DHS_mapped_path, TF_mapped_path, 
                  epilogos_mapped_path):
    #top 100 region inds and scores
    rec_inds = np.argsort(-score_arr)[:100]
    top_region_inds = np.append(roi, rec_inds)
    top_region_scores = np.append(1, score_arr[rec_inds])

    #most contributing metrics to the top 100 scoring regions
    top_region_contributions = np.append('ROI', contributing_cosine_arr[rec_inds])
    top_region_contributions = np.where(top_region_contributions=='0', 'DHS', top_region_contributions)
    top_region_contributions = np.where(top_region_contributions=='1', 'TF', top_region_contributions)
    top_region_contributions = np.where(top_region_contributions=='2', 'epilogos', top_region_contributions)
    top_region_contributions = np.where(top_region_contributions=='3', 'Coexpression', top_region_contributions)
    top_region_contributions = np.where(top_region_contributions=='4', 'Cooccurrence', top_region_contributions)
    top_region_contributions = np.where(top_region_contributions=='5', 'Database', top_region_contributions)
    top_region_contributions = np.where(top_region_contributions=='6', 'Experiment', top_region_contributions)
    top_region_contributions = np.where(top_region_contributions=='7', 'Fusion', top_region_contributions)
    top_region_contributions = np.where(top_region_contributions=='8', 'Homology', top_region_contributions)
    top_region_contributions = np.where(top_region_contributions=='9', 'Neighborhood', top_region_contributions)
    top_region_contributions = np.where(top_region_contributions=='10', 'Textmining', top_region_contributions)


    
    #combining scores and contributions
    score_contribution_comb = np.char.add(np.char.add(np.round(top_region_scores, 4).astype(str), 
                                                      np.array([' - ']*len(top_region_scores))), 
                                          top_region_contributions)


    #adding top comp to roiurls
    #component order
    f = open(DHS_mapped_path, 'r')
    DHS_component_order = np.array(f.readline().strip().split('\t'))[3:]
    f.close()
    f = open(TF_mapped_path, 'r')
    TF_component_order = np.array(f.readline().strip().split('\t'))[3:]
    f.close()
    f = open(epilogos_mapped_path, 'r')
    epilogos_component_order = np.array(f.readline().strip().split('\t'))[3:]
    f.close()
    #ooe arrays
    DHS_ooe_hdf5 = h5py.File(DHS_ooe_path, 'r')
    DHS_ooe = DHS_ooe_hdf5['ooe_arr']
    TF_ooe_hdf5 = h5py.File(TF_ooe_path, 'r')
    TF_ooe = TF_ooe_hdf5['ooe_arr']
    epilogos_ooe_hdf5 = h5py.File(epilogos_ooe_path, 'r')
    epilogos_ooe = epilogos_ooe_hdf5['ooe_arr']
    #finding ooe of roi
    DHS_roi_rank = rankdata(-DHS_ooe[top_region_inds[0]])
    TF_roi_rank = rankdata(-TF_ooe[top_region_inds[0]])
    epilogos_roi_mean_rank = rankdata(-epilogos_ooe[top_region_inds[0]][:18])
    epilogos_roi_max_rank = rankdata(-epilogos_ooe[top_region_inds[0]][18:])
    #ooe of top metric for each region and finding top factor for roi and each region
    metric_top_comps = ['']
    for ind in range(1, len(top_region_contributions)):
        metric = top_region_contributions[ind]
        if metric=='DHS':
            DHS_ind_rank = rankdata(-DHS_ooe[top_region_inds[ind]])
            top_comp = np.argsort((DHS_roi_rank+DHS_ind_rank)/2)[0]
            metric_top_comps.append(' (' + DHS_component_order[top_comp] + ')')
        elif metric=='TF':
            TF_ind_rank = rankdata(-TF_ooe[top_region_inds[ind]])
            top_comp = np.argsort((TF_roi_rank+TF_ind_rank)/2)[0]
            metric_top_comps.append(' (' + TF_component_order[top_comp] + ')')
        elif metric=='epilogos':
            epilogos_ind_mean_rank = rankdata(-epilogos_ooe[top_region_inds[ind]][:18])
            epilogos_ind_max_rank = rankdata(-epilogos_ooe[top_region_inds[ind]][18:])
            top_mean_comp = np.argsort((epilogos_roi_mean_rank + epilogos_ind_mean_rank)/2)[0]
            top_max_comp = np.argsort((epilogos_roi_max_rank + epilogos_ind_max_rank)/2)[0]
            metric_top_comps.append(' (' + epilogos_component_order[top_mean_comp] + '/' + 
                                    epilogos_component_order[top_max_comp] + ')')
        else:
            metric_top_comps.append(' ')

    DHS_ooe_hdf5.close()
    TF_ooe_hdf5.close()
    epilogos_ooe_hdf5.close()
    #add to roiurl data
    metric_top_comps = np.array(metric_top_comps)
    score_contribution_comb = np.char.add(score_contribution_comb, metric_top_comps)


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
    roi_url_name = 'region_comparison_recommended_regions_roi_'+str(roi)+'_scale_25kb_CT' + today.strftime("%Y%m%d") + '.bed'
    bed_path_name = os.path.join('/home/ctrader/public_html/', roi_url_name)
    np.savetxt(bed_path_name, top_regions, delimiter='\t', fmt='%s')
    link_path = 'https://epilogos.altius.org:3001/?application=viewer&sampleSet=vC&mode=single&genome=hg38&model=18&complexity=KL&group=all&chrLeft=chr1&chrRight=chr1&start=8000&stop=35000&roiURL=https%3A%2F%2Fencode%3Acollabor8%40resources.altius.org%2F~ctrader%2F' + roi_url_name    
    return link_path, bed_path_name

def load_hdf5_files(hdf5_path, roi):
    #take in a file path and load scaled scores and similar regions
    hdf5_file = h5py.File(hdf5_path, 'r')
    scaled_vals = hdf5_file['top_1000_scaled_vals'][roi]
    similar_regions = hdf5_file['top_1000_similar_regions'][roi]
    hdf5_file.close()
    return scaled_vals, similar_regions

def main(gz_file_path, genomic_location, DHS_hdf5_path, TF_hdf5_path, epilogos_hdf5_path, coex_hdf5_path, cooccur_hdf5_path, database_hdf5_path, exper_hdf5_path, 
	 fusion_hdf5_path, homology_hdf5_path, neigh_hdf5_path, textmin_hdf5_path, regions_path, DHS_ooe_path, TF_ooe_path, epilogos_ooe_path, 
         DHS_mapped_path, TF_mapped_path, epilogos_mapped_path):
    t = time.time()

    #find roi
    gz_tb = tabix.open(gz_file_path)
    roi = int(list(gz_tb.querys(genomic_location))[0][3])

    print('roi = ' + str(roi), flush=True)

    ###Load rank and ind arrays
    #DHS counts
    DHS_count_scaled_vals, DHS_count_similar_regions = load_hdf5_files(DHS_hdf5_path, roi)
    #TF motif counts
    TF_count_scaled_vals, TF_count_similar_regions = load_hdf5_files(TF_hdf5_path, roi)
    #epilogos score sum and max
    epilogos_score_scaled_vals, epilogos_score_similar_regions = load_hdf5_files(epilogos_hdf5_path, roi)
    #coexpression STRING scores
    coex_scaled_vals, coex_similar_regions = load_hdf5_files(coex_hdf5_path, roi)
    #cooccurence STRING scores
    cooccur_scaled_vals, cooccur_similar_regions = load_hdf5_files(cooccur_hdf5_path, roi)
    #database STRING scores
    database_scaled_vals, database_similar_regions = load_hdf5_files(database_hdf5_path, roi)
    #experiments STRING scores
    exper_scaled_vals, exper_similar_regions = load_hdf5_files(exper_hdf5_path, roi)
    #fusion STRING scores
    fusion_scaled_vals, fusion_similar_regions = load_hdf5_files(fusion_hdf5_path, roi)
    #homology STRING scores
    homology_scaled_vals, homology_similar_regions = load_hdf5_files(homology_hdf5_path, roi)
    #neighborhood STRING scores
    neigh_scaled_vals, neigh_similar_regions = load_hdf5_files(neigh_hdf5_path, roi)
    #textmining STRING scores
    textmin_scaled_vals, textmin_similar_regions = load_hdf5_files(textmin_hdf5_path, roi)

    ###Combine scores 
    score_arr, contributing_cosine_arr = combine_scores(DHS_count_scaled_vals, DHS_count_similar_regions, TF_count_scaled_vals, TF_count_similar_regions, epilogos_score_scaled_vals, 
                                                        epilogos_score_similar_regions, coex_scaled_vals, coex_similar_regions,	cooccur_scaled_vals, cooccur_similar_regions, 
							database_scaled_vals, database_similar_regions, exper_scaled_vals, exper_similar_regions, fusion_scaled_vals, fusion_similar_regions, 
							homology_scaled_vals, homology_similar_regions, neigh_scaled_vals, neigh_similar_regions, textmin_scaled_vals, textmin_similar_regions)

    ###Return regions
    bed_path_name = region_return(score_arr, contributing_cosine_arr, roi, regions_path, DHS_ooe_path, TF_ooe_path, epilogos_ooe_path, DHS_mapped_path, TF_mapped_path, 
                                  epilogos_mapped_path)
    print(bed_path_name)
    print(time.time()-t)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], 
         sys.argv[13], sys.argv[14], sys.argv[15], sys.argv[16], sys.argv[17], sys.argv[18], sys.argv[19], sys.argv[20])
