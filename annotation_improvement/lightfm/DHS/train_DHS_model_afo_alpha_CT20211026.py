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
    a_val = int(argv)/10000
    train_model(a_val)

def train_model(a_val):
    dhs = pd.read_table('/home/ctrader/AG/AugmentedGenomics/recommendersystems/lightfm/DHS/dat_bin_FDR01_hg38.txt.gz', header=None).T
    data_coo = sp.sparse.coo_matrix(dhs.values)
    
    #80/20 train test split to test precision values between train and test datasets
    (train, test) = cross_validation.random_train_test_split(data_coo, test_percentage=0.2, random_state=np.random.RandomState(3))
    data_80_20 = {'train': train, 'test': test}
    #creating model with L2 regularization set to a_val, k=16 components
    model = LightFM(no_components=16, loss='warp', item_alpha=a_val, user_alpha=a_val, learning_rate=0.01)
    #training model
    print('Training model for 15 epochs')
    model.fit(data_80_20['train'], epochs=15, num_threads=4)

    model_path = os.path.join('/home/ctrader/Hold/DHS/', 'DHS_model_e15_d16_reg'+str(a_val)+'_lr01_rs3.pickle')
    print('model path =', model_path)

    data_path = os.path.join('/home/ctrader/Hold/DHS/', 'DHS_data_e15_d16_reg'+str(a_val)+'_lr01_rs3.pickle')
    print('data path =', data_path)

    pickle_out_model = open(model_path, 'wb')
    pickle.dump(model, pickle_out_model)
    pickle_out_model.close()

    pickle_out_data = open(data_path, 'wb')
    pickle.dump(data_80_20, pickle_out_data)
    pickle_out_data.close()
    

if __name__ == "__main__":
    main(sys.argv[1])
