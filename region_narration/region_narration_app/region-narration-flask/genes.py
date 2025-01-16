from flask import Blueprint, request
from common import hilbert_pos_to_order
# import os
import time
import json
from constants import full_gene_path
from common import tabix_query, hilbert_pos_to_order

# Create a Blueprint
genes_bp = Blueprint('genes_bp', __name__)


"""
Collects genes overlapping each segment
"""
def get_genes_for_regions(regions):
    start = time.time()
    print(f"Collecting genes for {len(regions)} regions", flush=True)

    genes = []
    for r in regions:
        if "coordinates" in r:
            r_coordinates = r["coordinates"]
        else:
            r_order = int(r["order"])
            r_chrm = r["chromosome"]
            r_i = r["i"]
            r_i_start = hilbert_pos_to_order(r_i, r_order, 14)
            r_i_end = hilbert_pos_to_order(r_i + 1, r_order, 14)
            r_coordinates = f"{r_chrm}:{r_i_start}-{r_i_end}"

        # attach genes overlapping region
        r_genes = tabix_query(r_coordinates, full_gene_path)
        r_genes_parsed = [
            {
                'chromosome': str(g[0]), 
                'domain_start': int(g[1]),
                'domain_end': int(g[2]),
                'name': str(g[3]),
                'stable_ID': str(g[4]),
                'strand': int(g[5]),
                'start': int(g[6]),
                'end': int(g[7]),
                'description': str(g[10]),
                'in_gene': bool(
                    (
                        ((r_i_start <= int(g[6])) and (r_i_end > int(g[6]))) or 
                        ((r_i_start < int(g[7])) and (r_i_end >= int(g[7]))) or 
                        ((r_i_start > int(g[6])) and (r_i_end < int(g[7])))
                    )
                )
            } for g in r_genes
        ]
        genes.append(r_genes_parsed)

    print(f"Time to collect genes: {time.time() - start}", flush=True)
    return json.dumps({"genes": genes})


"""
Endpoint to collect genes for each segment
"""
@genes_bp.route("fetch_genes", methods=['POST'])
def fetch_genes():

    request_info = request.get_json()
    regions = request_info['regions'] if 'regions' in request_info else None # [{chromosome, order, i}, ...] or [{coordinates}, ...]

    # collect genes for each segment
    genes = get_genes_for_regions(regions)

    return genes