import pandas as pd
import numpy as np
from lightfm import LightFM
import pickle
import sys
import os


def precision_recall_setup(model, data):
    """
    Function to set up precision and recall calculations. Takes in model and data. Returns a Numpy array of DHS
    recommendations for each biosample and a Numpy array of known positive interatcions between each biosample and
    DHS. 
    """
    data_array = data.toarray()
    n_biosamples, n_DHSs = data.shape
    recs = []
    kp = []    
    for i in range(n_biosamples): #n_biosamples
        biosample_rec_scores = model.predict((i), np.arange(n_DHSs))
        recs.append(np.argsort(-biosample_rec_scores))
        kp.append(set(np.where(data_array[i] == 1)[0]))
    recs = np.array(recs)
    kp = np.array(kp)
    return recs, kp


def precision_recall_calc(data, recs, kp, k_vals=[1, 10, 100, 1000, 10000]):
    """
    Function to calculate precision recall values. Takes in a Numpy array of DHS recommendations for each 
    biosample, a Numpy array of known positive interatcions between each biosample and DHS, and a list of
    k values. Returns a list of the average precision values across all biosamples and a list of the 
    average recall values across all biosamples
    """
    precision_list = []
    recall_list = []
    n_biosamples, n_DHSs = data.shape
    for k in k_vals:
        top_k_recs = []
        for i in range(n_biosamples): #n_biosamples
            top_k_recs.append(set(recs[i, :k]))
        top_k_recs = np.array(top_k_recs)
        kp_topk = kp-top_k_recs
        topk_kp = top_k_recs-kp
        recs_len = np.repeat(k, n_biosamples) #n_biosamples
        fn = []
        fp = []
        for i in range(n_biosamples): #n_biosamples
            fn.append(len(kp_topk[i]))
            fp.append(len(topk_kp[i]))
        fn = np.array(fn)
        fp = np.array(fp)
        tp = recs_len - fp
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        precision_list.append(np.average(precision))
        recall_list.append(np.average(recall))
        print('precision / recall@%s = %s / %s' %(k, np.average(precision), np.average(recall)))
    return precision_list, recall_list

def precision_recall_full(model, data):
    """
    Full precision recall calculation function. Takes in model, data, and a list of k values. Returns a list of 
    the average precision values across all biosamples and a list of the average recall values across all biosamples.
    """
    print('Running precision and recall calculation setup')
    recs, kp = precision_recall_setup(model, data)
    print('Calculating precision and recall')
    precision, recall = precision_recall_calc(data, recs, kp)
    return precision, recall


def main(model_name, data_name, train_or_test,  k_vals):
    model_path = os.path.join('/home/ctrader/Hold/DHS/', model_name+'.pickle')
    pickle_in = open(model_path, 'rb')
    model=pickle.load(pickle_in)
    data_path = os.path.join('/home/ctrader/Hold/DHS/', data_name+'.pickle')
    pickle_in = open(data_path, 'rb')
    data=pickle.load(pickle_in)
    precision, recall = precision_recall_full(model, data[train_or_test])
    write_path = os.path.join('/home/ctrader/AG/AugmentedGenomics/recommendersystems/lightfm/DHS/precision_recall_results/', model_name+train_or_test+'precision_recall.txt')
    file = open(write_path, 'w')
    file.write("%s = %s \n%s = %s" %('Precision@'+str(k_vals), precision, 'Recall@'+str(k_vals), recall))
    file.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

