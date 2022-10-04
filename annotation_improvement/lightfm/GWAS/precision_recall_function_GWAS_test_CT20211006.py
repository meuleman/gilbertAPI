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
    pickle_in = open('/home/ctrader/Hold/GWAS/gwas_model_e30_d10_reg00001_CT20211006.pickle', 'rb')
    model = pickle.load(pickle_in)
    pickle_in = open('/home/ctrader/Hold/GWAS/gwas_data_e30_d10_reg00001_CT20211006.pickle', 'rb')
    data = pickle.load(pickle_in)
    func = partial(precision_recall, model, data['test'])
    with closing(Pool(5)) as pool:
        results = pool.map(func, [1, 10, 100, 1000, 10000])
    pool.join()
    precision = []
    recall = []
    for i in results:
        p, r = i
        precision.append(p)
        recall.append(r)
    file = open("GWAS_model_e30_d10_reg00001_test_precision_recall.txt", "w")
    file.write("%s = %s \n%s = %s" %("Test Precision@[1, 10, 100, 1000, 10000]", precision, "Test Recall@[1, 10, 100, 1000, 10000]", recall))
    file.close()
    
def precision_recall(model, data, k):
    #Takes in the model, data as a sparse coo matrix, and the values of k
    #for calculating precision and recall at k. Returns the calculated 
    #precision and recall for the model and data as a tuple.
    dense = data.todense()
    n_traits, n_snps = data.shape
    print('k = ' + str(k))
    num = []
    n_known_pos = []
    for trait in range(n_traits):
        if trait % 500 == 0:
            print(trait, '/', n_traits, 'traits done')
        scores = model.predict(trait, np.arange(n_snps)) 
        top_snps = [np.argsort(-scores)][0][0:k]
        known_pos = [i for i, snp_pos in enumerate(np.squeeze(np.asarray(dense[trait]))) if snp_pos > 0]
        top_in_known = len([i for i in known_pos if i in top_snps])
        num.append(top_in_known)
        known_pos_len = len(known_pos)
        if (known_pos_len == 0) and (top_in_known == 0):
            known_pos_len = 1
        n_known_pos.append(known_pos_len)
            
    num = np.array(num)
    precision = (num/k).mean()
    recall = (num/n_known_pos).mean()
    print('precision@%s = %s' %(k, precision))
    print('recall@%s = %s' %(k, recall))
    return precision, recall
    
if __name__ == "__main__":
    main()

