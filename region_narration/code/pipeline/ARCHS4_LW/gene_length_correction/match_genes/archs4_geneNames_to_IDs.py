import pandas as pd
import pickle
import csv

# Open ARCHS4 gene namesfile
#root = "/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW"
root = "/Users/lwu/PycharmProjects/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW"
pickle_in = open(root + "/Analysis/checkUnmatchedGenes_readCountSums/archs4_geneSymbolNames_list.pkl", 'rb')
archs4_geneNames_list = pickle.load(pickle_in)
pickle_in.close()

# Open the g:Profiler gene symbols to gene names converter
entire_data = pd.read_csv('gProfiler_archs4_geneSymbol_geneIDs.csv')
data = entire_data.iloc[:,:3] # ['initial_alias', 'converted_alias', 'name']; Don't need the 4th column ('Description')
data = data.to_numpy()

sv_data = []
unsaved = 0
for gene in data:
    gene_init = gene[0]
    gene_id = gene[1]
    gene_alias = gene[2]

    if gene_init == gene_alias:
        sv_gene = [gene_id, gene_init]
        sv_data.append(sv_gene)
    elif gene_init.startswith("ENSG"):
        sv_gene = [gene_init, gene_init] # Sometimes, gene_id is "None" when initial input was a gene ID
        sv_data.append(sv_gene)
    else:
        print(gene_init, "/t", gene_alias)
        unsaved += 1

print()
print("Data Size:", len(data))
print("Unsaved:", unsaved)
print("Saved:", len(sv_data))
print()
print(sv_data[0:5])

with open('archs4_geneSymbols_to_IDs.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(sv_data)


