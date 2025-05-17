import pandas as pd
import numpy as np
from lightfm import LightFM
import pickle
import sys
import os


def pr(model, test_data, train_data, k_vals):
    n_rows, n_cols = test_data.shape
    test_data_np = test_data.toarray()
    train_data_np = train_data.toarray()
    ordered_list = []
    test_count_list = []
    for row in range(n_rows):
        #known positives in train set
        train_kp = np.where(train_data_np[row]==1)[0]
        #every DHS index
        full_ind = np.arange(n_cols)
        #only the DHS indices that are not known positives in the train set
        test_only_ind = np.delete(full_ind, train_kp)
        #For each biosample, the DHS values for each DHS that does not exist as known positive in train set
        test_only_data = test_data_np[row, test_only_ind]
        #For each biosample, counting the number of true positives in test set only
        test_only_count = np.sum(test_only_data)
        test_count_list.append(test_only_count)
        #For each biosample, predicting scores for only dhs indices that are not known positives in train set
        scores = model.predict(row, test_only_ind)
        #ordering scores
        ordered = np.argsort(-scores)
        #ordering DHS values that were not known positives in train set 
        ordered_test_only = test_only_data[ordered]
        ordered_list.append(ordered_test_only)
    precision = []
    recall = []
    for k in k_vals:
        k_precision = []
        k_recall = []
        for row in range(n_rows):
            #for each biosample, calculates precision for the top k DHS recs
            count = np.sum(ordered_list[row][:k])
            k_precision.append(count/k) #this is how lightFM calculates precision
            k_recall.append(count/test_count_list[row])
        precision.append(k_precision)
        recall.append(k_recall)
    return precision, recall

def main(model_name, data_name, train_or_test):
    if train_or_test == 'test':
        other_train_test = 'train'
    else:
        other_train_test = 'test'
    model_path = os.path.join('/home/ctrader/Hold/DHS/', model_name+'.pickle')
    pickle_in = open(model_path, 'rb')
    model=pickle.load(pickle_in)
    data_path = os.path.join('/home/ctrader/Hold/DHS/', data_name+'.pickle')
    pickle_in = open(data_path, 'rb')
    data=pickle.load(pickle_in)
    n_biosamples, n_DHSs = data['test'].shape
    k_vals = list(range(1, 1000))
    k_vals.extend(list(np.arange(1000, n_DHSs, 10000)))
    precision_recall = pr(model, data[train_or_test], data[other_train_test], k_vals)
    write_path = os.path.join('/home/ctrader/AG/AugmentedGenomics/recommendersystems/lightfm/DHS/precision_recall_results/', model_name+'_'+train_or_test+'_precision_recall_CT20211028.pickle')
    pickle_out_pr = open(write_path, 'wb')
    pickle.dump(precision_recall, pickle_out_pr)
    pickle_out_model.close()
    #file = open(write_path, 'w')
    #file.write("%s = %s \n%s = %s" %('Precision@'+str(k_vals), precision, 'Recall@'+str(k_vals), recall))
    #file.close()
    
    
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
