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
    
    pickle_in = open('/home/ctrader/Hold/GWAS/pickled_subthreshold_data/top_100_SNP_recommendations_for_each_trait_all_traits_id_list_CT20220124.pickle', 'rb')
    top_100_recommended_snps_all_traits_id_list=pickle.load(pickle_in)

    find_top_100_SNPs(trait_values[trait_num:trait_num+5], trait_num, top_100_recommended_snps_all_traits_id_list)

def find_top_100_SNPs(tv, trait_num, top_100_recommended_snps_all_traits_id_list):
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
    for original_trait in tv:
        trait = original_trait.replace(')', ' ').replace('(', ' ').replace('-', ' ').replace('[', ' ').replace(']', ' ').replace('+', ' ').replace('_', ' ').replace(' ', '')
        if len(trait_codes_df[trait_codes_df['description']==(trait)])==1:
            phenotype_list.append(list(trait_codes_df[trait_codes_df['description']==(trait)]['phenotype'])[0])
        elif ((len(trait_codes_df[trait_codes_df['description']==(trait)])==2) & (len(trait_codes_df[(trait_codes_df['description']==(trait)) & (trait_codes_df['phenotype'].str.contains('irnt'))])==1)):
            phenotype_list.append(list(trait_codes_df[(trait_codes_df['description']==(trait)) & (trait_codes_df['phenotype'].str.contains('irnt'))]['phenotype'])[0])
        else:
            phenotype_list.append('SKIP_CT20220107')
            not_added.append(original_trait)

    #search through full data to find the p-values for each recommended SNP for each trait
    pvals_for_rec_snps = []
    snp_id_check = []
    for phenotype, trait, i in zip(phenotype_list, tv, range(len(phenotype_list))):
        i = i + trait_num
        if phenotype != 'SKIP_CT20220107':
            t = time.time()
            print(trait, phenotype)
            trait_data = []
            with gzip.open('/home/nasi4/proj/GWAS/UKBB_201807/phenotypes/' + phenotype + '.gwas.imputed_v3.both_sexes/' + phenotype + '.hg38.bed.gz','r') as f:
                for line in f:
                    trait_data.append(line.decode("utf-8").strip('\n').split('\t'))
            trait_df = pd.DataFrame(trait_data, columns=['Chromosome', 'Start', 'End', 'Reference SNP cluster ID', 'p-value'])
            trait_df['p-value'] = trait_df['p-value'].astype('float')

            rec_snp_ids_one_trait = top_100_recommended_snps_all_traits_id_list[(i*100):(i*100)+100]
        
            trait_df_snp_ind_based = trait_df.set_index('Reference SNP cluster ID')
            trait_df_snp_ind_based = trait_df_snp_ind_based.loc[trait_df_snp_ind_based.index.isin(rec_snp_ids_one_trait)]
            trait_df_snp_ind_based = trait_df_snp_ind_based.sort_values(['p-value'])
            trait_df_snp_ind_based = trait_df_snp_ind_based.groupby(['Reference SNP cluster ID']).first()
        
            if len(trait_df_snp_ind_based)==100:
                trait_df_snp_ind_based = trait_df_snp_ind_based.loc[rec_snp_ids_one_trait]
                pvals_for_rec_snps.extend(list(trait_df_snp_ind_based['p-value']))
                snp_id_check.extend(list(trait_df_snp_ind_based.index))
        
            else:
                trait_df_snp_ind_based['SNP ID Col 2'] = trait_df_snp_ind_based.index
                for snp_id in rec_snp_ids_one_trait:
                    snp_df = trait_df_snp_ind_based[trait_df_snp_ind_based['SNP ID Col 2']==snp_id]
                    snp_id_check.append(snp_id)
                    if len(snp_df)==1:
                        pvals_for_rec_snps.append(snp_df['p-value'])
                    else:
                        pvals_for_rec_snps.append(np.nan)
                        
            print(time.time()-t)
        else:
            pvals_for_rec_snps.extend(list(np.array(['SKIP_CT20220107']*100)))
            snp_id_check.extend(list(np.array(['SKIP_CT20220107']*100)))

    #pickle out the p-values and traits 
    pickle_out = open('/home/ctrader/Hold/GWAS/pickled_subthreshold_data/recommended_SNPs_for_each_trait_'+str(trait_num)+'_through_'+str(trait_num+4)+'_CT20220124.pickle', 'wb')
    pickle.dump(snp_id_check, pickle_out)
    pickle_out.close()

    pickle_out = open('/home/ctrader/Hold/GWAS/pickled_subthreshold_data/recommended_SNP_pvals_for_each_trait_'+str(trait_num)+'_through_'+str(trait_num+4)+'_CT20220124.pickle', 'wb')
    pickle.dump(pvals_for_rec_snps, pickle_out)
    pickle_out.close()

    pickle_out = open('/home/ctrader/Hold/GWAS/pickled_subthreshold_data/phenotype_list_for_recs_'+str(trait_num)+'_through_'+str(trait_num+4)+'_CT20220124.pickle', 'wb')
    pickle.dump(phenotype_list, pickle_out)
    pickle_out.close()

if __name__ == "__main__":
    main(sys.argv[1])
