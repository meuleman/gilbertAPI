import pandas as pd
import pickle

# Open ARCHS4 gene-symbol list
root = "/Users/lwu/PycharmProjects/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW"
file_path = "/Analysis/checkUnmatchedGenes_readCountSums/archs4_geneSymbolNames_list.pkl"
pickle_in = open(root + file_path, 'rb')
geneNames_list = pickle.load(pickle_in)
pickle_in.close()
print(len(geneNames_list))

# Open ARCHS4 gene names to ID file
df = pd.read_csv(root + "/gene_length_correction/match_genes/archs4_TLF_geneLengthData.csv") # Has header

# Extract values from the second column ('gene_name')
gene_name_column = df.iloc[:,1]
print(len(gene_name_column))
print(gene_name_column[0:5])

# Convert the column to a list
gene_name_list = gene_name_column.tolist()

matched_gene_idxs = []
unmatched_gene_idxs = []
for idx, gene_name in enumerate(geneNames_list):
    if gene_name.upper() in gene_name_list: # Here we are referring to the second column of the numpy array
        matched_gene_idxs.append(idx)
    else:
        unmatched_gene_idxs.append(idx)
        #print(gene_name)

# for idx, gene_name in enumerate(gene_name_list):
#     if gene_name.upper() in geneNames_list: # Here we are referring to the second column of the numpy array
#         matched_gene_idxs.append(idx)
#     else:
#         unmatched_gene_idxs.append(idx)
#         #print(gene_name)


print()
print("Complete")
print("Matched:", len(matched_gene_idxs))
#print("ie:", matched_gene_idxs[0:5])
print("Unmatched:", len(unmatched_gene_idxs))
#
# print(matched_gene_idxs[0], matched_gene_idxs[-1])

pickle_out = open("archs4_matched_geneIdx_LW20230120.pkl", 'wb')
pickle.dump(matched_gene_idxs, pickle_out)
pickle_out.close()

pickle_out = open("archs4_unmatched_geneIdx_LW20230120.pkl", 'wb')
pickle.dump(unmatched_gene_idxs, pickle_out)
pickle_out.close()
