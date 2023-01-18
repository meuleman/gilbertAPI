import json
import pickle
import h5py

# Open ARCHS4 v2.1 Data
infile = "/home/ajoshi/archs4/data/archs4_gene_human_v2.1.h5"

# - Use read_direct method to improve read performance
with h5py.File(infile, 'r') as hdf5_file:
    data = hdf5_file['data']
    expression = data['expression'][:]

# Open the JSON file
root = "/home/lwu/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/sbatch_files"
with open('/results/allSelected_uniqueSampleIndices_list.json') as json_file:
    # Load the JSON data from the file
    sampleIdxs_list = json.load(json_file)

# Extract rows from the matrix using the indexes in the list
extracted_matrix = expression[:, sampleIdxs_list]

pickle_out = open("root + /results/samples_extracted_matrix.pkl", 'wb')
pickle.dump(extracted_matrix, pickle_out)
pickle_out.close()