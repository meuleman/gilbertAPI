from itertools import islice
from collections import Counter
import numpy as np
import h5py
import re
import json

# Open ARCHS4 v2.1 Data
infile = "/home/ajoshi/archs4/data/archs4_gene_human_v2.1.h5"

# - Use read_direct method to improve read performance
with h5py.File(infile, 'r') as hdf5_file:
    data = hdf5_file['data']
    expression = data['expression'][:]

    meta = hdf5_file['meta']
    symbols = meta['genes']['gene_symbol'][:]
    source_name_ch1 = meta['samples']['source_name_ch1'][:]

# Decode source names using list comprehension
source_decoded = [x.decode('UTF-8') for x in source_name_ch1]

# Clean source names using regular expressions
source_cleaned = [re.sub(r'[-/\'\\ _]', '', x) for x in source_decoded]

# List of all categories (keywords)
category_list = [  # Cardiovascular
    "Heart", "Atrium", "Cardiac Muscle Fiber", "Pericardium", "Valve", "Ventricle",
    "HUVEC",  # Endothelial

    # Connective Tissue
    "Adipose", "Adipocyte",
    "Myofibroblast", "Osteoblast", "Stromal Cell", "MG63", "K562",  # Bone
    "Bone Marrow",
    "Chondrocyte",
    "HT1080", "NHDF",  # Connective

    # Digestive
    "Esophagus",
    "Colon", "Colonic Muscosa", "Intestinal Epithelial Cell", "Small Intestine",
    "Liver", "Hepatic Stellate Cell", "Hepatocyte", "Primary Liver Cell", "HEPG2",  # Liver
    "Alpha Cell", "Beta Cell", "Pancreatic Islet", "PANC1",  # Pancreas
    "Gastric Epithelial Cell",
    "HCT116", "HT29", "RKO",  # Connective, Colon

    # Immune
    "Granulocyte", "Neutrophil",
    "Blymphocyte", "Plasma Cell", "Plasmacytoid Dendritic Cell", "Tlymphocyte",
    "JURKAT", "MT4", "NALM6", "REH", "RPMI8226",  # Lymphoid
    "Alveolar", "Dendritic", "Kupffer", "Macrophage", "Macroglia",
    "THP1", "U937",  # Myeloid
    "JEG3",  # Placenta
    "Spleen",
    "Thymus", "Thymocyte",

    # Integumentary
    "Skin", "Basal Cell", "Fibroblast", "Hair Follicle", "Keratinocyte", "Melanocyte",
    "A431", "HACAT", "HNSCC",

    # Muscular
    "Skeletal Muscle", "Myoblast",
    "Airway Smooth Muscle",  # "Myofibroblast", [Repeat]
    "Respiratory Smooth Muscle", "Vascular Smooth Muscle",
    "LHCNM2",  # Muscle

    # Nervous
    "Astrocyte", "Cerebellum", "Cerebral Cortex", "Hypothalamus", "Medula",
    "Midbrain", "Neuron", "Oligodendrocyte", "Pons", "Spinal Cord", "Thalamus",
    "Retina",
    "Motor Neuron", "Sensory Neuron",
    "Brain",  # Example in the ARCHS4 Website [1]
    "DAOY", "IMR32", "SKNSH",  # Brain

    # Respiratory
    # "Alveolar", [Repeat]
    "Lung", "Lung Epithelial",
    "A549", "H1299", "IMR90", "MRC5", "NCIH460", "NHBE",  # Lung
    "Trachea",

    # Urogenital
    "Breast", "Mammary Gland", "MCF7", "MCF10", "MDAMB231", "MDAMB468", "SKBR3", "T47D",
    "HELA",  # Cervix
    "Kidney", "Podocyte", "Renal Cortex",
    "293F", "FLPIN TREX 293", "HEK293",  # Kidney
    "C42", "DU145", "LNCAP", "PC3",  # Prostate

    "Ovary", "Granulosa", "Oocyte",
    "ES2", "SKOV3",  # Ovary
    "Testis",
    "GSE81547",  # Example in the ARCHS4 GitHub
    "GSE30017", "GSM2679484"]  # Examples in the ARCHS4 Website [2]

# Create a dictionary that maps categories to their sample indices
used_indices = set()
category_sampleIdx_dict = {}
for category in category_list:
    category_sampleIdx_dict[category] = []
    for idx, string in enumerate(source_cleaned):
        if (category.replace(" ", "").lower() in string.lower()) and (idx not in used_indices):
            used_indices.add(idx)
            category_sampleIdx_dict[category].append(idx)

# Get the total sample size collected for each keyword
category_sampleSizes_dict = {key: len(category_sampleIdx_dict[key]) for key in category_sampleIdx_dict}

# Flatten the list values into a single list of integers
allSelected_uniqueSampleIndices_list = list(set([val for vals in category_sampleIdx_dict.values() for val in vals]))

root = "~/work/NarrativeGenomics/region_narration/code/pipeline/ARCHS4_LW/python_files"
with open("/results/category_sampleIdx_dict.json", "w") as f:
    json.dump(category_sampleIdx_dict, f)

with open("/results/category_sampleSizes_dict.json", "w") as f:
    json.dump(category_sampleSizes_dict, f)

with open("/results/allSelected_uniqueSampleIndices_list.json", "w") as f:
    json.dump(allSelected_uniqueSampleIndices_list, f)