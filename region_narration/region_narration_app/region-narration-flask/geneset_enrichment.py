import json
from flask import Blueprint, request
import numpy as np
from constants import gsXgene, genes, genesets
from common import perform_geneset_hypergeometric_test, generate_geneset_return_data


# Create a Blueprint
geneset_enrichment_bp = Blueprint('geneset_enrichment_bp', __name__)


# Helper functions
def get_geneset_genes(geneset):
    """
    Collects the genes in a geneset
    
    Args:
        geneset (str): The name of the geneset
        
    Returns:
        list: List of genes in the geneset
        
    Raises:
        ValueError: If geneset is not found
    """
    try:
        geneset_genes = genes[gsXgene[genesets == geneset][0]]
        return geneset_genes.tolist()

    except IndexError:
        raise ValueError(f"Geneset '{geneset}' not found")


# Route handlers
@geneset_enrichment_bp.route("geneset_enrichment", methods=['POST'])
def geneset_enrichment():
    """
    Endpoint for calculating enrichment for genesets given a list of genes
    
    Request:
        genes (list): List of gene names
        threshold (float, optional): P-value threshold, default 0.05
        membership (bool, optional): Whether to return membership, default False
        bonferroni (bool, optional): Whether to apply Bonferroni correction, default False
        
    Returns:
        JSON: List of enriched genesets with statistics
    """
    try:
        request_info = request.get_json()
        if not request_info or "genes" not in request_info:
            return json.dumps({"error": "Missing required parameter: genes"}), 400
        
        genes_for_test = request_info["genes"]
        threshold = request_info.get("threshold", 0.05)
        membership = str(request_info.get("membership", "False")).lower() == "true"
        bonferroni = str(request_info.get("bonferroni", "False")).lower() == "true"

        # test enrichment of genes (if membership, p is bool array)
        p = perform_geneset_hypergeometric_test(
            genes_for_test, gsXgene, genes, 
            genesets, membership=membership, 
            bonferroni=bonferroni
        )

        # generate return data
        return_list = generate_geneset_return_data(
            p, genesets, gsXgene, 
            genes, threshold=threshold
        )
        return json.dumps(return_list)

    except Exception as e:
        return json.dumps({"error": str(e)}), 500


@geneset_enrichment_bp.route("geneset_to_genes", methods=['POST'])
def geneset_to_genes():
    """
    Endpoint for inverse geneset search, where a geneset is provided and a list of genes is returned
    
    Request:
        geneset (str): The name of the geneset to query
        
    Returns:
        JSON: List of genes in the specified geneset
    """
    try:
        request_info = request.get_json()
        if not request_info or "geneset" not in request_info:
            return json.dumps({"error": "Missing required parameter: geneset"}), 400
        
        # get genes in geneset
        geneset = request_info['geneset']
        geneset_genes = get_geneset_genes(geneset)
        return json.dumps(geneset_genes)
    
    except ValueError as e:
        return json.dumps({"error": str(e)}), 404
    
    except Exception as e:
        return json.dumps({"error": str(e)}), 500