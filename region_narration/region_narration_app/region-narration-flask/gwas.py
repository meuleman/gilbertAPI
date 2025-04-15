from flask import Blueprint, request
import numpy as np
import requests
import json
import h5py
import os
import time
from constants import gwas_dir as datadir

# Create a Blueprint
gwas_bp = Blueprint('gwas_bp', __name__)

# preload gwas data for faster access
with h5py.File(os.path.join(datadir, "ukbb.94_traits.0.5_pip.h5"), 'r') as f:
    datasets = list(f.keys())
    gwas_dict = {f[key][:] for key in datasets}

"""
Fetches the traits associated with any genomic position
"""
def fetch_gwas(chromosome, position):
    dataset_name = f"{chromosome}:{position}"
    # with h5py.File(os.path.join(datadir, "ukbb.94_traits.0.5_pip.h5"), 'r') as f:
    try:
        data = gwas_dict[dataset_name][:]
        # scores
        scores = data[1].tolist()
        # traits
        trait_indices = data[0].astype(int)
        unique_trait_indices = np.unique(trait_indices)
        # collect trait names
        unique_trait_names = gwas_dict['traits'][unique_trait_indices].astype(str)
        index2trait = {i: t for i, t in zip(unique_trait_indices, unique_trait_names)}
        trait_names = list(map(index2trait.get, trait_indices))
    except:
        trait_names = []
        scores = []
    return trait_names, scores


"""
Endpoint for collecting the GWAS associations for any set of order 14 single bp positions
"""
@gwas_bp.route("gwas_for_positions", methods=['POST'])
def gwas_for_positions():
    request_info = request.get_json()
    positions = request_info['positions']  # [{chromosome, index}, ...]

    return_object = []
    for position in positions:
        # collect GWAS associations for positions one at a time
        trait_names, scores = fetch_gwas(position['chromosome'], position['index'])
        return_object.append(
            {
                "chromosome": position['chromosome'], 
                "index": position['index'], 
                "trait_names": trait_names, 
                "scores": scores
            }
        )
    print(return_object)

    return json.dumps(return_object)