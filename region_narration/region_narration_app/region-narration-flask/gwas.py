from flask import Blueprint, request
import numpy as np
import requests
import json
import h5py
import os
import time

# Create a Blueprint
gwas_bp = Blueprint('gwas_bp', __name__)

datadir="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/gwas/20240823"

"""
Fetches the traits associated with any genomic position
"""
def fetch_gwas(chromosome, position):
    dataset_name = f"{chromosome}:{position}"
    with h5py.File(os.path.join(datadir, "sig_associations.v1.0.1e-8.h5"), 'r') as f:
        try:
            data = f[dataset_name][:]
            # scores
            scores = data[1].tolist()
            # traits
            trait_indices = data[0].astype(int)
            unique_trait_indices = np.unique(trait_indices)
            # collect trait names
            unique_trait_names = f['traits'][unique_trait_indices]
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

    return json.dumps(return_object)