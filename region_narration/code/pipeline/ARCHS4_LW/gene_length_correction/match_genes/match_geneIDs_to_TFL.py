import pandas as pd
import numpy as np
import csv

# Open and read the gene-length file
with open("transcriptLengthPerGene_WM20221126.txt", "r") as f:
    lines = f.readlines()[:]

genes_data = [line.strip().split("\t") for line in lines]
wm_genes_symbol, wm_genes_id, wm_gene_length = zip(*genes_data)
wm_dict = {gene_symbol: gene_length for gene_symbol, _, gene_length in genes_data}

# Open ARCHS4 gene names to ID file
data = pd.read_csv('archs4_geneSymbols_to_IDs.csv', header=None)
data = data.to_numpy()

sv_data = []
unsaved_ = []
unsaved = 0
diff = 0
for gene in data:
    gene_id = gene[0]
    gene_name = gene[1]

    if gene_id in wm_genes_id:
        index = wm_genes_id.index(gene_id)
        #print(index)
        sv_gene = [gene_id, gene_name, wm_gene_length[index]]
        sv_data.append(sv_gene)

        if gene_name != wm_genes_symbol[index]:
            #print(gene_name, "|", wm_genes_symbol[index])
            diff += 1

    else:
        #print(gene_id, "/t", gene_name)
        unsaved_.append(gene_name)
        unsaved += 1

print()
print("Complete")
print("Data Size:", len(data))
print("Unsaved:", unsaved)
print("Saved:", len(sv_data))
print("Different Names:", diff)
print()
print(sv_data[0:5])

with open('archs4_TLF_geneLengthData.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(sv_data)