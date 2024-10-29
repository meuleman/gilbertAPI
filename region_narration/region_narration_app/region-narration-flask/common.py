import h5py
import numpy as np
import time
import os


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


"""
Function to preload all indices to be used for filter previews up through the max order
"""
def preload_preview_indices(index_path, max_order=12):
    preload_time = time.time()
    preloaded_indices = {}
    
    with h5py.File(index_path, 'r') as f:

        # find all the datasets less than or equal to the max order
        datasets = list(f.keys())
        datasets = [d for d in datasets if int(d.split('.')[0]) <= max_order]
        
        # preload indices for each dataset
        for d in datasets:

            # create chromosome independent key
            order, chrm, dataset_name, factor_index = d.split('.')
            key = f"{order}.{dataset_name}.{factor_index}"

            # all factor instances across chromosomes are saved together in a list
            if key not in preloaded_indices:
                preloaded_indices[key] = []
            
            # add to factorXorder list
            preloaded_indices[key].append({
                "chromosome": chrm,
                "order": int(order), 
                "index": int(factor_index), 
                "dataset_name": dataset_name, 
                "indices": f[d][:, 0].astype("int32"), 
                "scores": f[d][:, 1].astype("float32")
            })

    # save the order through which indices have been preloaded
    preloaded_indices['max_order'] = max_order
            
    print(f"Time to preload indices for filter previews: {time.time() - preload_time}", flush=True)
    return preloaded_indices


"""
Converts an array of hilbert positions from one order to another
"""
def hilbert_pos_to_order(pos, frm, to):
    return np.floor(pos / pow(4, frm - to)).astype(int)