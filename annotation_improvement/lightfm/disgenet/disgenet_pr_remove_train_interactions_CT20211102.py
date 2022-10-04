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
        train_kp = np.where(train_data_np[row]>0)[0]
        #every gene index
        full_ind = np.arange(n_cols)
        #only the gene indices that are not known positives in the train set
        test_only_ind = np.delete(full_ind, train_kp)
        #For each disease, the gene association values for each gene that does not exist as known positive in train set
        test_only_data = test_data_np[row, test_only_ind]
        #For each disease, counting the number of true positives in test set only
        test_only_count = np.count_nonzero(test_only_data > 0)
        test_count_list.append(test_only_count)
        #For each disease, predicting scores for only gene indices that are not known positives in train set
        scores = model.predict(row, test_only_ind)
        #ordering scores
        ordered = np.argsort(-scores)
        #ordering gene values that were not known positives in train set 
        ordered_test_only = test_only_data[ordered]
        ordered_list.append(ordered_test_only)
    precision = []
    recall = []
    for k in k_vals:
        k_precision = []
        k_recall = []
        for row in range(n_rows):
            #for each disease, calculates precision for the top k gene recs
            count = np.count_nonzero(ordered_list[row][:k] > 0)
            k_precision.append(count/k) #this is how lightFM calculates precision
            if count == 0 and test_count_list[row] == 0:
                k_recall.append(0)
            else:
                k_recall.append(count/test_count_list[row])
        precision.append(k_precision)
        recall.append(k_recall)
    return precision, recall

def main(model_name, data_name, train_or_test):
    if train_or_test == 'test':
        other_train_test = 'train'
    else:
        other_train_test = 'test'
    model_path = os.path.join('/home/ctrader/Hold/Dis/', model_name+'.pickle')
    pickle_in = open(model_path, 'rb')
    model=pickle.load(pickle_in)
    data_path = os.path.join('/home/ctrader/Hold/Dis/', data_name+'.pickle')
    pickle_in = open(data_path, 'rb')
    data=pickle.load(pickle_in)
    n_diseases, n_genes = data['test'].shape
    k_vals = [1, 10, 100, 1000, n_genes]
    precision_recall = pr(model, data[train_or_test], data[other_train_test], k_vals)
    write_path = os.path.join('/home/ctrader/AG/AugmentedGenomics/recommendersystems/lightfm/disgenet/precision_recall_results/', model_name+'_'+train_or_test+'_precision_recall_CT20211028.pickle')
    pickle_out_pr = open(write_path, 'wb')
    pickle.dump(precision_recall, pickle_out_pr)
    pickle_out_pr.close()
    #file = open(write_path, 'w')
    #file.write("%s = %s \n%s = %s" %('Precision@'+str(k_vals), precision, 'Recall@'+str(k_vals), recall))
    #file.close()
    
    
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
