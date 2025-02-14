from flask import Blueprint, request
import numpy as np
import json
import h5py
import os
import time
from constants import gsXgene, genes, genesets
from common import perform_geneset_hypergeometric_test, generate_geneset_return_data


# Create a Blueprint
geneset_enrichment_bp = Blueprint('geneset_enrichment_bp', __name__)


def get_geneset_genes(geneset):
    """
    Collects the genes in a geneset
    """
    geneset_genes = genes[gsXgene[genesets == geneset][0]]
    return geneset_genes.tolist()


@geneset_enrichment_bp.route("geneset_enrichment", methods=['POST'])
def geneset_enrichment():
    """
    Endpoint for calculating enrichment for genesets given a list of genes
    """
    request_info = request.get_json()
    genes_for_test = request_info["genes"]  # ["gene1", "gene2", ...]
    threshold = request_info.get("threshold", 0.05)  # default threshold is 0.05
    membership = str(request_info.get("membership", "False")).lower() == "true"  # default membership is False, convert to boolean
    bonferroni = str(request_info.get("bonferroni", "False")).lower() == "true"  # default bonferroni is False, convert to boolean

    # test enrichment of genes (if membership, p is bool array)
    p = perform_geneset_hypergeometric_test(genes_for_test, gsXgene, genes, genesets, membership=membership, bonferroni=bonferroni)

    # generate return data
    return_list = generate_geneset_return_data(p, genesets, gsXgene, genes, threshold=threshold)

    return return_list


@geneset_enrichment_bp.route("geneset_to_genes", methods=['POST'])
def geneset_to_genes():
    """
    Endpoint for inverse geneset search, where a geneset is provided and a list of genes is returned
    """
    request_info = request.get_json()
    geneset = request_info['geneset']

    # get genes in geneset
    geneset_genes = get_geneset_genes(geneset)

    return json.dumps(geneset_genes)