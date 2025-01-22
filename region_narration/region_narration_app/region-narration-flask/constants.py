import numpy as np
import os
from common import preload_preview_indices

### DATA DIRECTORIES ###
filtering_path="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/data_filtering/20241101/indices.scores.native_resolution.h5"
order_independent_filtering_path="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/order_independent_filtering/20241205/overlap_merged.h5"
csn_dir="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/precomputed_csn/20241205"
partial_csn_dir="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/precomputed_partial_csn/20241219"
gwas_dir="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/gwas/20241021"
genedir="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/genes/20240906/"
full_gene_path = os.path.join(
    genedir, 
    "mart_export.protein_coding.basal_plus_extension.5kb_upstream.1kb_downstream.1000kb_distal.CT20240906.bed.bgz"
)
# columns for full_gene_path
# [
# 'Chromosome', 'Left distal', 'Right distal', 'Gene name', 
# 'Gene stable ID', 'Strand', 'Gene start (bp)', 'Gene end (bp)', 
# 'Basal upstream', 'Basal downstream', 'Gene description'
# ]

# Path to binary matrix file
geneset_dir="/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/geneset_enrichment/20240918"
geneset_bin_path = os.path.join(geneset_dir, "c5.go.bp.v2024.1.bin.h5")
########################

### CONSTANTS ###
CHR_SIZES_UCSC = {
    "chr1": 248956422,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr2": 242193529,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chrX": 156040895,
    "chrY": 57227415
}


ORDER_SEGMENT_COUNTS = {
    4: 2957,
    5: 11792,
    6: 47136,
    7: 188505,
    8: 753983,
    9: 3015901,
    10: 12063568,
    11: 48254229,
    12: 193016875,
    13: 772067468,
    14: 3088269856
}


chromosomes = np.append([f"chr{c}" for c in range(1, 23)], ['chrX', 'chrY'])
orders = np.arange(4, 15)
precalculated_region_set_enrichment_orders = [4, 5]
#################

### DATA ###

# preload all indices for filter previews
PRELOADED_INDICES = preload_preview_indices(filtering_path)
############


####PUBMED SUMMARY####
MODEL_ID = "nomic-ai/nomic-embed-text-v1.5"
LANCEDB_URI = "/data/narrativeGenomics/NarrativeGenomics/region_narration/data/gilbert/pubmed_summary/20250122/lancedb/jcr_mini"