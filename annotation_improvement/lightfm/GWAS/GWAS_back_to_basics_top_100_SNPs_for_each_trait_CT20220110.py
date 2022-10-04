import sys
import time
import pandas as pd
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import pickle
import gzip
import textwrap
import os
from datetime import date

def main(argv1):
    trait_num = int(argv1)
    pickle_in = open('/home/ctrader/Hold/GWAS/GWAS_trait_values_CT20220110.pickle', 'rb')
    trait_values=pickle.load(pickle_in)
    find_top_100_SNPs(trait_values[trait_num:trait_num+5], trait_num)

def find_top_100_SNPs(trait_values, trait_num):
    #load traits and phenotypes
    trait_codes = []
    with gzip.open('/home/nasi4/proj/GWAS/UKBB_201807/info/phenotypes.both_sexes.tsv.gz','r') as f:        
        for line in f:
            trait_codes.append(line.decode("utf-8").strip('\n').split('\t'))

    trait_codes_df = pd.DataFrame(trait_codes[1:], columns=trait_codes[0])

    trait_codes_df['description'] = trait_codes_df['description'].str.replace(')', ' ')
    trait_codes_df['description'] = trait_codes_df['description'].str.replace('(', ' ')
    trait_codes_df['description'] = trait_codes_df['description'].str.replace('-', ' ')
    trait_codes_df['description'] = trait_codes_df['description'].str.replace('[', ' ')
    trait_codes_df['description'] = trait_codes_df['description'].str.replace(']', ' ')
    trait_codes_df['description'] = trait_codes_df['description'].str.replace('+', ' ')
    trait_codes_df['description'] = trait_codes_df['description'].str.replace('_', ' ')
    trait_codes_df['description'] = trait_codes_df['description'].str.replace(' ', '')

    #find the phenotype for each trait
    phenotype_list = []
    not_added = []
    for original_trait in trait_values:
        trait = original_trait.replace(')', ' ').replace('(', ' ').replace('-', ' ').replace('[', ' ').replace(']', ' ').replace('+', ' ').replace('_', ' ').replace(' ', '')
        if len(trait_codes_df[trait_codes_df['description']==(trait)])==1:
            phenotype_list.append(list(trait_codes_df[trait_codes_df['description']==(trait)]['phenotype'])[0])
        elif ((len(trait_codes_df[trait_codes_df['description']==(trait)])==2) & (len(trait_codes_df[(trait_codes_df['description']==(trait)) & (trait_codes_df['phenotype'].str.contains('irnt'))])==1)):
            phenotype_list.append(list(trait_codes_df[(trait_codes_df['description']==(trait)) & (trait_codes_df['phenotype'].str.contains('irnt'))]['phenotype'])[0])
        else:
            phenotype_list.append('SKIP_CT20220107')
            not_added.append(original_trait)

    #search through full data to find the top 100 SNPs for each trait
    top_100_SNPs_for_each_trait = []
    for phenotype, trait in zip(phenotype_list, trait_values):
        if phenotype != 'SKIP_CT20220107':
            t = time.time()
            print(trait, phenotype)
            trait_data = []
            with gzip.open('/home/nasi4/proj/GWAS/UKBB_201807/phenotypes/' + phenotype + '.gwas.imputed_v3.both_sexes/' + phenotype + '.hg38.bed.gz','r') as f:        
                for line in f:        
                    trait_data.append(line.decode("utf-8").strip('\n').split('\t'))
            trait_df = pd.DataFrame(trait_data, columns=['Chromosome', 'Start', 'End', 'Reference SNP cluster ID', 'p-value'])
            trait_df['p-value'] = trait_df['p-value'].astype('float')
            trait_df = trait_df[trait_df['p-value']!=0]
            trait_df = trait_df.sort_values(['p-value']).reset_index()
            top_100_SNPs_for_each_trait.append(trait_df.iloc[:100])
            print(time.time()-t)
        else:
            top_100_SNPs_for_each_trait.append('SKIP_CT20220107')

    #pickle out the full data, the phenotypes used, and the traits not included
    pickle_out = open('/home/ctrader/Hold/GWAS/pickled_subthreshold_data/top_100_SNPs_for_each_trait_'+str(trait_num)+'_through_'+str(trait_num+4)+'_CT20220110.pickle', 'wb')
    pickle.dump(top_100_SNPs_for_each_trait, pickle_out)
    pickle_out.close()

    pickle_out = open('/home/ctrader/Hold/GWAS/pickled_subthreshold_data/phenotype_list_'+str(trait_num)+'_through_'+str(trait_num+4)+'_CT20220110.pickle', 'wb')
    pickle.dump(phenotype_list, pickle_out)
    pickle_out.close()

    pickle_out = open('/home/ctrader/Hold/GWAS/pickled_subthreshold_data/not_added_'+str(trait_num)+'_through_'+str(trait_num+4)+'_CT20220110.pickle', 'wb')
    pickle.dump(not_added, pickle_out)
    pickle_out.close()

if __name__ == "__main__":
    main(sys.argv[1])
