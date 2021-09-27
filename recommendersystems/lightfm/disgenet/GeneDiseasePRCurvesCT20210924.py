#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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


def main():
    pickle_in = open('gd_model100.pickle', 'rb')
    model=pickle.load(pickle_in)
    
    pickle_in = open('gd_data100.pickle', 'rb')
    data=pickle.load(pickle_in)
    
    # measure precision/recall @ k for k=1:50
    test_precision = []
    test_recall = []
    train_precision = []
    train_recall = []

    for i in np.arange(1, 51):
        test_precision.append(precision_at_k(model, data['test'], k=i).mean())
        test_recall.append(recall_at_k(model, data['test'], k=i).mean())
        train_precision.append(precision_at_k(model, data['train'], k=i).mean())
        train_recall.append(recall_at_k(model, data['train'], k=i).mean())
        print('k=' + str(i) + ' Done')
    
    #Plot precision-recall curves for train and test data
    plt.plot(test_recall, test_precision)
    plt.plot(train_recall, train_precision)
    plt.xlabel('Average Recall@k')
    plt.ylabel('Average Precision@k')
    plt.legend(['Test Set', 'Train Set'])
    plt.title('Precision-Recall Curves for Disease vs Gene LightFM model')
    plt.savefig('Precision_Recall_Curves_for_Disease_vs_Gene_LightFM_model.jpeg')

if __name__ == "__main__":
    main()

