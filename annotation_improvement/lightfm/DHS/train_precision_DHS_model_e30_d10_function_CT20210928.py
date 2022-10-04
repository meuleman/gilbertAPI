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
    pickle_in = open('DHS_model30.pickle', 'rb')
    model=pickle.load(pickle_in)
    pickle_in = open('DHS_data30.pickle', 'rb')
    data=pickle.load(pickle_in)
    func = partial(precision_fun, model, data['train'])
    with closing(Pool(5)) as pool:
        results = pool.map(func, [1, 10, 100, 1000, 10000])
    pool.join()
    file = open("DHS_model_e30_d10_train_precision.txt", "w")
    file.write("%s = %s" %("Train Precision@[1, 10, 100, 1000, 10000]", results))
    file.close()
    
def precision_fun(model, data, k_val):
    #precision function
    return precision_at_k(model, data, k=k_val).mean()
    
if __name__ == "__main__":
    main()


# In[ ]:




