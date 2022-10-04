#!/usr/bin/env python
# coding: utf-8

# In[2]:

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
    model=pickle.load(pickle_in)
    pickle_in = open('/home/ctrader/Hold/DHS/DHS_data_e15_d16_reg0001.pickle', 'rb')
    data=pickle.load(pickle_in)
    func = partial(precision_fun, model, data['test'], data['train'])
    with closing(Pool(5)) as pool:
        results = pool.map(func, [1, 10, 100, 1000, 10000])
    pool.join()
    pickle_out_results = open('/home/ctrader/Hold/DHS/precision_CT20211028.pickle', 'wb')
    pickle.dump(results, pickle_out_results)
    pickle_out_results.close()
    #file = open("test_only_precision.txt", "w")
    #file.write("%s = %s" %("Test Precision@[1, 10, 100, 1000, 10000]", results))
    #file.close()
    
def precision_fun(model, test_data, train_data, k_val):
    #precision function
    return precision_at_k(model, test_data, train_interactions=train_data, k=k_val, check_intersections=False)
    
if __name__ == "__main__":
    main()


# In[ ]:




