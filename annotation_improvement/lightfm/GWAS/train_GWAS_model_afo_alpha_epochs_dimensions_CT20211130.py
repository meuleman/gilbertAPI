import sys
import os
import pandas as pd
import scipy as sp
import numpy as np
from lightfm import cross_validation
from lightfm import LightFM
import matplotlib.pyplot as plt
import pickle
from datetime import date

def main(argv1, argv2, argv3):
    a_val = int(argv1)/100000
    num_epochs = int(argv2)
    num_dimensions = int(argv3)
    model=train_model(a_val, num_epochs, num_dimensions)

def train_model(a_val, num_epochs, num_dimensions):
    #load in pickled data split into train (80%), tune (10%), and test (10%)
    pickle_in = open('/home/ctrader/Hold/GWAS/pickled_models_afo_alpha_epochs_dimensions/gwas_data_for_training_CT20211130.pickle', 'rb')
    gwas_data=pickle.load(pickle_in)

    #creating model with L2 regularization set to a_val, k=num_dimensions components
    model = LightFM(no_components=num_dimensions, loss='warp', item_alpha=a_val, user_alpha=a_val, learning_rate=0.01)
    #training model
    print('Training model for' + str(num_epochs) + 'epochs')
    model.fit(gwas_data['train'], epochs=num_epochs, num_threads=4)


    today = date.today()
    model_path = os.path.join('/home/ctrader/Hold/GWAS/pickled_models_afo_alpha_epochs_dimensions/', 'GWAS_model_e'+str(num_epochs)+'_d'+str(num_dimensions)+'_reg'+str(a_val)+'_lr01_CT'+today.strftime("%Y%m%d")+'.pickle')
    print('model path =', model_path)

    pickle_out_model = open(model_path, 'wb')
    pickle.dump(model, pickle_out_model)
    pickle_out_model.close()
    
    return model

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
