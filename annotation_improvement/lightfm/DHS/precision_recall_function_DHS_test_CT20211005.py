#!/usr/bin/env python
# coding: utf-8

# In[ ]:

from functools import partial
import pandas as pd
import scipy as sp
import numpy as np
from lightfm import cross_validation
from lightfm.evaluation import precision_at_k
from lightfm.evaluation import recall_at_k
from lightfm.evaluation import auc_score
from lightfm import LightFM
import matplotlib.pyplot as plt
import pickle
from multiprocessing import Pool
from contextlib import closing


def main():
    #read in model and data
    pickle_in = open('/home/ctrader/Hold/DHS/DHS_model_e15_d16_reg0001.pickle', 'rb')
    model = pickle.load(pickle_in)
    pickle_in = open('/home/ctrader/Hold/DHS/DHS_data_e15_d16_reg0001.pickle', 'rb')
    data = pickle.load(pickle_in)
    func = partial(precision_recall, model, data['test'])
    with closing(Pool(8)) as pool:
        results = pool.map(func, [1, 10, 100, 1000, 10000, 100000, 1000000, 3591898])
    pool.join()
    precision = []
    recall = []
    for i in results:
        p, r = i
        precision.append(p)
        recall.append(r)
    file = open("/precision_recall_results/DHS_model_e15_d16_reg0001_test_precision_recall.txt", "w")
    file.write("%s = %s \n%s = %s" %("Test Precision@[1, 10, 100, 1000, 10000, 100000, 1000000, 3591898]", precision, "Test Recall@[1, 10, 100, 1000, 10000, 100000, 1000000, 3591898]", recall))
    file.close()
    
def precision_recall(model, data, k):
    #Takes in the model, data as a sparse coo matrix, and the values of k
    #for calculating precision and recall at k. Returns the calculated 
    #precision and recall for the model and data as a tuple of lists.
    dense = data.todense()
    n_biosamples, n_dhs = data.shape
    print('k = ' + str(k))
    num = []
    n_known_pos = []
    for biosample in range(733):
        if biosample % 100 == 0:
            print(biosample, '/', n_biosamples, 'biosamples done')
        scores = model.predict(biosample, np.arange(n_dhs)) 
        top_DHSs = [np.argsort(-scores)][0][0:k]
        known_pos = [i for i, dhs_pos in enumerate(np.squeeze(np.asarray(dense[biosample]))) if dhs_pos == 1]
        top_in_known = len([i for i in known_pos if i in top_DHSs])
        num.append(top_in_known)
        n_known_pos.append(len(known_pos))
            
    num = np.array(num)
    precision = (num/k).mean()
    recall = (num/n_known_pos).mean()
    print('precision@%s = %s' %(k, precision))
    print('recall@%s = %s' %(k, recall))
    return precision, recall
    
if __name__ == "__main__":
    main()

