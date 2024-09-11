import h5py
import numpy as np


"""
Fetches the order 14 positional indices for a given order, chromosome, factor index, 
and dataset name. Returns an empty array if the indices are not found.
"""
def fetch_indices_scores(
    order, 
    chromosome, 
    factor_index, 
    dataset_name, 
    file_path
):
    filter_key = f"{order}.{chromosome}.{dataset_name}.{factor_index}"
    with h5py.File(file_path, 'r') as f:
        try:
            indices_scores = f[filter_key][:]
        except:
            indices_scores = np.array([], dtype="int32")
    return indices_scores