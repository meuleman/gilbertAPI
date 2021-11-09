import sys
import os
import pandas as pd
import scipy as sp
import numpy as np
from lightfm import cross_validation
from lightfm import LightFM
import matplotlib.pyplot as plt
import pickle

def main(argv):
    a_val = int(argv)/100000
    train_model(a_val)

def train_model(a_val):
    #read in data
    tsv_file = open("curated_gene_disease_associations.tsv")
    df = pd.read_csv(tsv_file, delimiter="\t")
    gd_df = df[['geneId', 'diseaseId', 'score']]
    #finding unique genes and diseases
    geneId_values = np.unique(gd_df['geneId'].values)
    diseaseId_values = np.unique(gd_df['diseaseId'].values)
    #Creating sparse matrix from csv formatted data
    gd_matrix = pd.DataFrame(index = diseaseId_values, columns = geneId_values)
    for disease in diseaseId_values:
        disease_df = gd_df[gd_df['diseaseId']==disease].set_index('geneId')
        gd_matrix.loc[disease] = disease_df['score']
    gd_np = gd_matrix.to_numpy()
    gd_np = np.array(gd_np, dtype=float)
    gd_np[np.isnan(gd_np)] = 0
    gd_coo = sp.sparse.coo_matrix(gd_np)


    #80/20 train test split to test precision values between train and test datasets
    (train, test) = cross_validation.random_train_test_split(gd_coo, test_percentage=0.2, random_state=None)
    data_80_20 = {'train': train, 'test': test, 'gene_ids' : geneId_values}
    #creating model with L2 regularization set to a_val, k=16 components
    model = LightFM(no_components=10, loss='warp', item_alpha=a_val, user_alpha=a_val, learning_rate=0.05)
    #training model
    print('Training model for 100 epochs')
    model.fit(data_80_20['train'], epochs=100, num_threads=4)

    model_path = os.path.join('/home/ctrader/Hold/Dis/', 'disgenet_model_e100_d10_reg'+str(a_val)+'_lr05.pickle')
    print('model path =', model_path)

    data_path = os.path.join('/home/ctrader/Hold/Dis/', 'disgenet_data_e100_d10_reg'+str(a_val)+'_lr05.pickle')
    print('data path =', data_path)

    pickle_out_model = open(model_path, 'wb')
    pickle.dump(model, pickle_out_model)
    pickle_out_model.close()

    pickle_out_data = open(data_path, 'wb')
    pickle.dump(data_80_20, pickle_out_data)
    pickle_out_data.close()
    

if __name__ == "__main__":
    main(sys.argv[1])
