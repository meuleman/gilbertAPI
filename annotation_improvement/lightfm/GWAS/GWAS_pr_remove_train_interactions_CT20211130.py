import pandas as pd
import numpy as np
from lightfm import LightFM
import pickle
import sys
import os
from datetime import date


def pr(model, tune_data, train_data, test_data, k_vals):
    n_rows, n_cols = tune_data.shape
    tune_data_np = tune_data.toarray()
    #tune_data_np = np.ceil(tune_data_np).astype(int)
    #tune_data_np = np.where(tune_data_np > 0, 1, 0)
    train_test_data = train_data + test_data
    train_test_data_np = train_test_data.toarray()
    #train_test_data_np = np.where(train_test_data_np > 0, 1, 0)
    #train_test_data_np = np.ceil(train_test_data_np).astype(int)
    ordered_list = []
    tune_count_list = []
    for row in range(n_rows):
        #known positives in train and test sets
        train_test_kp = np.where(train_test_data_np[row]>0)[0]
        #every SNP index
        full_ind = np.arange(n_cols)
        #only the SNP indices that are not known positives in the train and test sets
        tune_only_ind = np.delete(full_ind, train_test_kp)
        #For each trait, the SNP values for each SNP that does not exist as known positive in train or test set
        tune_only_data = tune_data_np[row, tune_only_ind]
        #For each trait, counting the number of true positives in tune set only
        tune_only_count = np.count_nonzero(tune_only_data > 0)
        tune_count_list.append(tune_only_count)
        #For each trait, predicting scores for only SNP indices that are not known positives in train set
        scores = model.predict(row, tune_only_ind)
        #ordering scores
        ordered = np.argsort(-scores)
        #ordering SNP values that were not known positives in train or test sets 
        ordered_tune_only = tune_only_data[ordered]
        ordered_list.append(ordered_tune_only)
    precision = []
    recall = []
    for k in k_vals:
        k_precision = []
        k_recall = []
        for row in range(n_rows):
            #for each trait, calculates precision for the top k SNP recs
            count = np.count_nonzero(ordered_list[row][:k] > 0)
            k_precision.append(count/k) #this is how lightFM calculates precision
            if (tune_count_list[row]==0) and (count==0):
                k_recall.append(0)
            else:
                k_recall.append(count/tune_count_list[row])
        precision.append(k_precision)
        recall.append(k_recall)
    return precision, recall

def main(model_name, data_name, train_tune_test):
    if train_tune_test == 'test':
        to_remove_1 = 'train'
        to_remove_2 = 'tune'
    elif train_tune_test == 'train':
        to_remove_1 = 'test'
        to_remove_2 = 'tune'
    else:
        to_remove_1 = 'train'
        to_remove_2 = 'test'
    model_path = os.path.join('/home/ctrader/Hold/GWAS/pickled_models_afo_alpha_epochs_dimensions/', model_name+'.pickle')
    pickle_in = open(model_path, 'rb')
    model=pickle.load(pickle_in)
    data_path = os.path.join('/home/ctrader/Hold/GWAS/pickled_models_afo_alpha_epochs_dimensions/', data_name+'.pickle')
    pickle_in = open(data_path, 'rb')
    data=pickle.load(pickle_in)
    n_traits, n_SNPs = data['tune'].shape
    k_vals = list(range(1, 1000))
    k_vals.extend(list(np.arange(1000, n_SNPs, 10000)))
    precision_recall = pr(model, data[train_tune_test], data[to_remove_1], data[to_remove_2], k_vals)
    today = date.today()
    write_path = os.path.join('/home/ctrader/Hold/GWAS/precision_recall_results_afo_alpha_epochs_dimensions/', model_name+'_'+train_tune_test+'_precision_recall_CT'+today.strftime("%Y%m%d")+'.pickle')
    pickle_out_pr = open(write_path, 'wb')
    pickle.dump(precision_recall, pickle_out_pr)
    pickle_out_pr.close()
    
    
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
