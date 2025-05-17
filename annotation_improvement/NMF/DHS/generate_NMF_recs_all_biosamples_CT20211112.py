import sys
import os
import pandas as pd
import scipy as sp
import numpy as np
import gzip
import matplotlib.pyplot as plt
import pickle
from multiprocessing import Pool
from contextlib import closing
from functools import partial

def main():
    dhs = pd.read_table('/home/ctrader/Hold/DHS/dat_bin_FDR01_hg38.txt.gz', header=None).T
    fn = gzip.GzipFile('/home/ctrader/Hold/DHS/2018-06-08NC16_NNDSVD_Basis.npy.gz', "r")
    decomp_Basis = np.load(fn)
    fn = gzip.GzipFile('/home/ctrader/Hold/DHS/2018-06-08NC16_NNDSVD_Mixture.npy.gz', "r")
    decomp_Mixture = np.load(fn)
    dhs_array = dhs.to_numpy()
    tsv_file = open("/home/ctrader/Hold/DHS/DHS_Index_and_Vocabulary_metadata.tsv")
    df = pd.read_csv(tsv_file, delimiter="\t")
    biosample_df = df[['library order', 'Biosample name']]
    DHS_metadata = pd.read_table('/home/ctrader/Hold/DHS/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz', header=None)

    col_names = DHS_metadata.iloc[0]
    DHS_metadata = DHS_metadata[1:]
    DHS_metadata.columns = col_names
    DHS_metadata = DHS_metadata.reset_index()
    DHS_metadata[['identifier']] = DHS_metadata[['identifier']].astype(str)
    DHS_metadata = DHS_metadata.drop(columns=['index'])    
    #reordering dataframe to match order of presence/absence matrix and DHS metadata dataframe
    f_eye_signal = pd.read_table('/home/ctrader/lightfm_data/out_ENCFF967FUD.txt', names=['name', 'size', 'covered', 'sum', 'mean0', 'mean'])
    f_eye_signal[['name']] = f_eye_signal[['name']].astype(str)
    f_eye_signal = f_eye_signal.set_index('name')
    id_order = DHS_metadata['identifier'].to_list()
    f_eye_signal = f_eye_signal.reindex(id_order)
    f_eye_signal = f_eye_signal.reset_index()



    func = partial(sample_recommendation, f_eye_signal, 0, decomp_Mixture, decomp_Basis, DHS_metadata, biosample_df, dhs, dhs_array, df)
    with closing(Pool(5)) as pool:
        rec_info_all = pool.map(func, list(range(733)))
    pool.join()
    pickle_out = open('/home/ctrader/AG/AugmentedGenomics/recommendersystems/NMF/DHS/pickled_data/rec_info_all_CT20211112.pickle', 'wb')
    pickle.dump(rec_info_all, pickle_out)
    pickle_out.close()



def signal_values(dhs, signal_data):
    identifier_series = pd.Series(signal_data['name'])
    identifier = identifier_series[dhs]
    #Average biosample signal across each recommended DHS
    mean_signal_sample_series = pd.Series(signal_data['mean'])
    mean_signal_sample = mean_signal_sample_series[dhs]
    #Average signal of all biosamples at each recommended DHS
    mean_signal_DHS_series = pd.Series(DHS_metadata['mean_signal'])
    mean_signal_DHS = mean_signal_DHS_series[dhs]
    return list(mean_signal_sample), list(mean_signal_DHS), list(identifier)


def sample_recommendation(signal_data, print_num, DHS_e, biosample_e, DHS_metadata, biosample_df, dhs, dhs_array, df, biosamples):
    
    #number of biosamples and DHSs in dataset
    n_biosamples, n_biosample_e = biosample_e.shape
    n_DHS_e, n_DHSs = DHS_e.shape
    return_list = []
    metadata = DHS_metadata.reset_index()
    comp_ser = metadata['component']
    chrom_ser = metadata['seqname']
    start_ser = metadata['start']
    end_ser = metadata['end']
    numsamples_ser = metadata['numsamples']
    identifier_ser = metadata['identifier']
    
    #generate DHS recommendations for each biosample
    for biosample_name in biosamples:
        if isinstance(biosample_name, str):
            biosample = list(biosample_df[biosample_df['Biosample name']==biosample_name]['library order'])[0] - 1
        else:
            biosample = biosample_name
        
        #DHSs known to relate to the biosample
        #indices
        known_positives = [i for i, dhs_pos in enumerate(list(dhs.T[biosample])) if dhs_pos == 1]
        
        #filter out known positives from what will be recommended
        not_known_ind = list(set(np.arange(n_DHSs)) - set(known_positives))
        
        nki_series = pd.Series(not_known_ind)
        
        #DHSs model predictions that relate to biosample
        filtered_scores = np.dot(decomp_Basis[biosample, :], decomp_Mixture[:, np.array(not_known_ind)])
        

        #rank from most related to least
        #not DHS indices
        filtered_top_DHSs = [np.argsort(-filtered_scores)][0]
        top_rec_DHSs = nki_series[filtered_top_DHSs]
        top_filtered_scores = filtered_scores[filtered_top_DHSs]
        
        comp_rec = list(comp_ser[top_rec_DHSs])
        chrom_rec = list(chrom_ser[top_rec_DHSs])
        start_rec = list(start_ser[top_rec_DHSs])
        end_rec = list(end_ser[top_rec_DHSs])
        numsamples_rec = list(numsamples_ser[top_rec_DHSs])
        identifier_rec = list(identifier_ser[top_rec_DHSs])
        
        count = 0
        total_count = 0
        return_biosample = []
        return_biosample.append(comp_rec)
        return_biosample.append(chrom_rec)
        return_biosample.append(start_rec)
        return_biosample.append(end_rec)
        return_biosample.append(numsamples_rec)
        return_biosample.append(list(top_rec_DHSs))
        return_biosample.append(list(top_filtered_scores))
        return_list.append(return_biosample)
        
        #print results
        print('Biosample: %s (%s)' %(biosample_name, str(biosample)))
        if print_num > 0:
            print('   Number of Known DHSs: %s' %(len(known_positives)))
            print()
            print('   Recommended DHS Locations: ') 
            print('   Location | Component | Number of Samples | Biosample Signal/Average Signal | Left +DHS Signal | Right +DHS Signal')
        
        
        
            dhs_ids = top_rec_DHSs.to_numpy()
            #dhs to the right of recommended dhs
            dhs_right = dhs_ids[:]+1
            #dhs to the left of recommended dhs
            dhs_left = dhs_ids[:]-1
        
            #removing dhss that are "neighbored" across chromosomes
            chrm_list = list(np.unique(list(DHS_metadata['seqname'])))
            rows_per_chrm = [0]
            for chrm in chrm_list:
                max_row = max(rows_per_chrm)
                rows, col = DHS_metadata[DHS_metadata['seqname']==chrm].shape
                rows_per_chrm.append(rows + max_row)
            rows_per_chrm = np.array(rows_per_chrm)
        
            dhs_right = np.delete(dhs_right, np.isin(dhs_right, rows_per_chrm))
            dhs_left = np.delete(dhs_left, np.isin(dhs_left, rows_per_chrm-1))
        
            #determining presence/absense of dhss to right and left of recommended dhs
            right = dhs_array[biosample, dhs_right] 
            left = dhs_array[biosample, dhs_left]
        
            #signal data
            sample_rec_dhs_signal, rec_dhs_signal, rec_identifier = signal_values(list(top_rec_DHSs), signal_data)
            sample_right_dhs_signal, right_dhs_signal, right_identifier = signal_values(list(dhs_right), signal_data)
            sample_left_dhs_signal, left_dhs_signal, left_identifier = signal_values(list(dhs_left), signal_data)
#         print(rec_identifier[:5])
#         print(identifier_rec[:5])
        
        
            for i in range(print_num):
#             print(identifier_rec[i])
#             print(rec_identifier[:5])
#             print(list(rec_identifier)[i])
                chrom = chrom_rec[i]
                start = start_rec[i]
                end = end_rec[i]
                component = comp_rec[i]
                num_samples = numsamples_rec[i]
                sample_rec_dhs_signal_i = round(list(sample_rec_dhs_signal)[i], 3)
                rec_dhs_signal_i = round(float(list(rec_dhs_signal)[i]), 3)
                dhs_id = list(top_rec_DHSs)[i]
                right_id = list(np.where(dhs_right == dhs_id+1))[0] #+1
                if list(right_id) and right[right_id]==1:
#                 print(int(right_id[0]))
                    sample_right_dhs_signal_i = round(sample_right_dhs_signal[int(right_id[0])], 3)
                    right_dhs_signal_i = round(right_dhs_signal[int(right_id[0])], 3)
                    right_print = str(round(sample_right_dhs_signal_i, 3)) + '/' + str(round(right_dhs_signal_i, 3))
                else:
                    right_print = 'N/A'
            
                left_id = list(np.where(dhs_left == dhs_id-1))[0]
                if list(left_id) and left[left_id]==1:
#                 print(int(left_id[0]))
                    sample_left_dhs_signal_i = round(sample_left_dhs_signal[int(left_id[0])], 3)
                    left_dhs_signal_i = round(left_dhs_signal[int(left_id[0])], 3)
                    left_print = str(round(sample_left_dhs_signal_i, 3)) + '/' + str(round(left_dhs_signal_i, 3))
                else:
                    left_print = 'N/A'
                
                
                print('      %s:%s-%s | %s | %s | %s/%s | %s | %s'
                      %(chrom, start, end, component, num_samples, sample_rec_dhs_signal_i, rec_dhs_signal_i,
                       left_print, right_print))

    return return_list



if __name__ == "__main__":
    main()
