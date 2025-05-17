import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
from datetime import date
import re
from tqdm import tqdm
import time
from pybedtools import BedTool
import random
import math
import pickle
import seaborn as sns
from scipy.stats import percentileofscore, pearsonr
from scipy.spatial.distance import jaccard, euclidean
from multiprocessing import Pool
from contextlib import closing
from functools import partial
import operator

today = date.today()



today = date.today()

gene_DHS_enrichment = pd.read_pickle('../data/gene-dhs-enrichment-dropna.pkl')

pearson_scores = []
    
def combined_func(gene_DHS_enrichment, gene):
    specific_gene_dhs_enrichment = gene_DHS_enrichment.loc[gene]
    
    my_pearsonr = partial(pearsonr, specific_gene_dhs_enrichment)

    gene_DHS_enrichment[['pearson_corr', 'p-val']] = pd.DataFrame(gene_DHS_enrichment.apply(my_pearsonr, axis=1).tolist(),
                                                                  index=gene_DHS_enrichment.index)
    
    return gene_DHS_enrichment['pearson_corr'].tolist()
    


func = partial(combined_func, gene_DHS_enrichment)
with closing(Pool(10)) as pool:
    results = pool.map(func, list(gene_DHS_enrichment.index.tolist()))
pool.join()

pearson_scores.extend(results) 
pearson_scores = np.array(pearson_scores)            
                   
pkl_file = open('../output/gene-dhs-distr.pkl','wb')
pickle.dump(pearson_scores, pkl_file)

fn = '../output/gene-dhs-distr.npy'
np.save(fn, pearson_scores)

distr_plot = sns.displot(pearson_scores, height=5, aspect=3)
distr_plot.set(title='Gene-DHS Enrichment Distribution')
plt.savefig('../output/plots/ST_' + today.strftime("%Y-%m-%d") + '_gene-DHS-enrichment-distribution.png')



# with open('../output/DHS_distribution.txt', 'w'):
    