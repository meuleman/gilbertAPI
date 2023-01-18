import h5py
import numpy as np

# # Open the HDF5 file
# with h5py.File('matrix.hdf5', 'r') as f:
#     # Read the matrix
#     matrix = f['matrix'][:]

matrix = np.array([[1,2,3,4,5,6,7],
                  [8,9,10,11,12,13,14],
                  [0,0,0,0,0,0,0]])

# List of unique integers
index_list = [0]

# Extract rows from the matrix using the indexes in the list
extracted_matrix = matrix[:, index_list]

print(extracted_matrix)
