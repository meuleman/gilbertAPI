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
    train_model(int(argv))

def train_model(epochs_val):
    dhs = pd.read_table('/home/ctrader/AG/AugmentedGenomics/recommendersystems/lightfm/DHS/dat_bin_FDR01_hg38.txt.gz', header=None).T
    data_coo = sp.sparse.coo_matrix(dhs.values)
    
    #80/20 train test split to test precision values between train and test datasets
    (train, test) = cross_validation.random_train_test_split(data_coo, test_percentage=0.2, random_state=None)
    data_80_20 = {'train': train, 'test': test}
    #creating model with L2 regularization set to 0.0001, k=16 components
    model = LightFM(no_components=16, loss='warp', item_alpha=0.0001, user_alpha=0.0001)
    #training model
    print('Training model for %s epochs' %(epochs_val))
    model.fit(data_80_20['train'], epochs=epochs_val, num_threads=4)

    model_path = os.path.join('/home/ctrader/Hold/DHS/', 'DHS_model_e'+str(epochs_val)+'_d16_reg0001.pickle')
    print('model path =', model_path)

    data_path = os.path.join('/home/ctrader/Hold/DHS/', 'DHS_data_e'+str(epochs_val)+'_d16_reg0001.pickle')
    print('data path =', data_path)

    pickle_out_model = open(model_path, 'wb')
    pickle.dump(model, pickle_out_model)
    pickle_out_model.close()

    pickle_out_data = open(data_path, 'wb')
    pickle.dump(data_80_20, pickle_out_data)
    pickle_out_data.close()
    

if __name__ == "__main__":
    main(sys.argv[1])
