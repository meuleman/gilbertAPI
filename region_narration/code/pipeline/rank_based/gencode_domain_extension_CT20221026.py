import numpy as np
import os
import sys
import re
from functools import partial
from multiprocessing import Pool
from contextlib import closing

def load_file(in_path):
    f = open(in_path, 'r')
    gencode_file = f.readlines()
    f.close()
    return gencode_file

def parse_data(gencode_file):
    data = []
    for line in gencode_file[5:]:
        line = re.split('\t|;| ', line.strip())
        #gene coordinates
        gene_chrm = line[0]
        gene_start = line[3]
        gene_end = line[4]
        gene_strand = line[6]
        gene_type = line[2]

        #gene name
        gene_name_index = line.index('gene_name')
        gene_name = line[gene_name_index + 1].strip("\"")

        #gene id (ENSG)
        gene_id_index = line.index('gene_id')
        gene_id = line[gene_id_index + 1].strip("\"")
        #protein id (ENSP)
        if 'protein_id' in line:
            protein_id_index = line.index('protein_id')
            protein_id = line[protein_id_index + 1].strip("\"")
        else:
            protein_id = 'NA'
        #promoter start and end
        if gene_strand == '+':
            promoter_start = gene_start
            promoter_end = str(int(gene_start) + 1)
            basal_start = str(int(gene_start) - 5000)
            basal_end = str(int(gene_start) + 1000)
            domain_start = str(int(gene_start) - 1000000)
            domain_end = str(int(gene_end) + 1000000)

        elif gene_strand == '-':
            promoter_start = str(int(gene_end) - 1)
            promoter_end = gene_end
            basal_start = str(int(gene_end) - 1000)
            basal_end = str(int(gene_end) + 5000)
            domain_start = str(int(gene_start) - 1000000)
            domain_end = str(int(gene_end) + 1000000)

        data_row = [gene_chrm, gene_start, gene_end, gene_name, gene_type, gene_strand, 
                    promoter_start, promoter_end, basal_start, basal_end, domain_start, 
                    domain_end, gene_id, protein_id]
        data.append(data_row)

    data = np.array(data)

    header = ['chrm', 'start', 'end', 'name', 'type', 'strand', 'Promoter start (bp)', 
              'Promoter end (bp)', 'Basal start (bp)', 'Basal end (bp)', 'Domain start (bp)', 
              'Domain end (bp)', 'gene id', 'protein id']
    return data, header

def create_mask_dicts(chrm, tss, domain_m, names):
    #creating chromosome masked data
    chrm_dom_dict = {}
    chrm_tss_dict = {}
    chrm_name_dict = {}

    for chromosome in np.unique(chrm):
        chrm_mask = chrm==chromosome
        chrm_doms = domain_m[chrm_mask, :]
        chrm_tss = tss[chrm_mask]
        chrm_names = names[chrm_mask]

        chrm_dom_dict[chromosome] = chrm_doms
        chrm_tss_dict[chromosome] = chrm_tss
        chrm_name_dict[chromosome] = chrm_names
    return chrm_dom_dict, chrm_tss_dict, chrm_name_dict

def limit_domains(chrm, tss, names, domain_m, chrm_dom_dict, chrm_tss_dict, chrm_name_dict, gene):
    gene_chrm = chrm[gene]
    gene_tss = tss[gene]
    gene_name = names[gene]
    gene_dom = domain_m[gene, :]
    
    #filtering out wrong chromosomes
    chrm_doms = chrm_dom_dict[gene_chrm]
    chrm_tss = chrm_tss_dict[gene_chrm]
    chrm_names = chrm_name_dict[gene_chrm]
    
    #filtering out other annotations of the same gene
    name_mask = chrm_names!=gene_name
    chrm_name_doms = chrm_doms[name_mask, :]
    chrm_name_tss = chrm_tss[name_mask]
    
    #find all the genes to the left and right
    left_genes = chrm_name_tss<gene_tss
    right_genes = chrm_name_tss>gene_tss
    left_doms = chrm_name_doms[left_genes, :]
    right_doms = chrm_name_doms[right_genes, :]

    #if there are any genes, find the closest basal region
    if len(left_doms) > 0:
        max_left_basal_end = np.max(left_doms[:, 1])
    else:
        max_left_basal_end = gene_dom[2]
    if len(right_doms) > 0:
        min_right_basal_start = np.min(right_doms[:, 0])
    else:
        min_right_basal_start = gene_dom[3]

    #setting the domain for each gene
    gene_dom[2] = np.min([np.max([max_left_basal_end, gene_dom[2]]), gene_dom[0]])
    gene_dom[3] = np.max([np.min([min_right_basal_start, gene_dom[3]]), gene_dom[1]])

    return gene_dom

def main(in_path, out_path, num_cpus):
    print('Loading...', flush=True)
    gencode_file = load_file(in_path)
    print('Parsing...', flush=True)
    data, header = parse_data(gencode_file)
    
    chrm = data[:, 0]
    tss = data[:, 6].astype(int)
    domain_m = data[:, 8:12].astype(int)
    names = data[:, 3]
    #limiting domain to real chromosome positions
    domain_m = np.where(domain_m < 1, 1, domain_m)

    print('Masks...', flush=True)
    chrm_dom_dict, chrm_tss_dict, chrm_name_dict = create_mask_dicts(chrm, tss, domain_m, names)

    print('Limiting domains...', flush=True)
    func = partial(limit_domains, chrm, tss, names, domain_m, chrm_dom_dict, chrm_tss_dict, chrm_name_dict)
    with closing(Pool(num_cpus)) as pool:
        new_domains = pool.map(func, np.arange(len(domain_m)))
    pool.join()
    new_domains = np.array(new_domains)

    #changing domain values
    data[:, 8:12] = new_domains.astype(str)

    print('Saving...', flush=True)
    np.savetxt(out_path, data, delimiter='\t', fmt='%s', header='\t'.join(header))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
