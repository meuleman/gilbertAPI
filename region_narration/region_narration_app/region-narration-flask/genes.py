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
def get_genes_for_segments(segments):
    start = time.time()
    print(f"Collecting genes for {len(segments)} regions", flush=True)

    genes = []
    for s in segments:
        s_order = int(s["order"])
        s_chrm = s["chromosome"]
        s_i = s["i"]
        s_i_start = hilbert_pos_to_order(s_i, s_order, 14)
        s_i_end = hilbert_pos_to_order(s_i + 1, s_order, 14)

        # attach genes overlapping region
        s_genes = tabix_query(f"{s_chrm}:{s_i_start}-{s_i_end}", full_gene_path)
        s_genes_parsed = [
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
                        ((s_i_start <= int(g[6])) and (s_i_end > int(g[6]))) or 
                        ((s_i_start < int(g[7])) and (s_i_end >= int(g[7]))) or 
                        ((s_i_start > int(g[6])) and (s_i_end < int(g[7])))
                    )
                )
            } for g in s_genes
        ]
        genes.append(s_genes_parsed)

    print(f"Time to collect genes: {time.time() - start}", flush=True)
    return json.dumps({"genes": genes})


"""
Endpoint to collect genes for each segment
"""
@genes_bp.route("fetch_genes", methods=['POST'])
def fetch_genes():

    request_info = request.get_json()
    regions = request_info['regions'] if 'regions' in request_info else None # [{chromosome, order, i}, ...]

    # collect genes for each segment
    genes = get_genes_for_segments(regions)

    return genes