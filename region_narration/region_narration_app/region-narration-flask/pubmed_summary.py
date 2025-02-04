from flask import Blueprint, request
# import numpy as np
import json
# import h5py
import os
import time
from sentence_transformers import SentenceTransformer
import lancedb
import outlines
from outlines import models, generate
import torch
import sys
# import pyarrow as pa
# from logger import ParquetRequestLogger
import uuid
# import datetime
from constants import MODEL_ID, LANCEDB_URI


# Create a Blueprint
pubmed_summary_bp = Blueprint('pubmed_summary_bp', __name__)


"""
Loads model
"""
def start_engine():
    device = torch.device("cpu")
    print("🥶 cold starting inference", flush=True)
    start = time.monotonic_ns()
    model = SentenceTransformer(MODEL_ID, trust_remote_code=True, device=device)
    model.eval()
    # print(f"CUDA memory allocated: {torch.cuda.memory_allocated() / 1e6} MB")
    duration_s = (time.monotonic_ns() - start) / 1e9
    print(f"🏎️ engine started in {duration_s:.0f}s", flush=True)
    db = lancedb.connect(LANCEDB_URI)
    table = db.open_table("jcr-mini")
    default_instructions = """"""
    examples = [
        {
            "query": "EWSR1/FLI1 MOTIF enrichment @ 16kbp; Stromal B DHS enrichment @ 1Mbp; Atrial fibrillation GWAS occurrence @ 1bp; Musculoskeletal DHS enrichment @ 64kbp; Cardiac DHS enrichment @ 256kbp; Quiescent/Low CS occurrence @ 256bp; NTMT2 GENE; GORAB GENE; N TERMINAL PROTEIN AMINO ACID MODIFICATION GO; EPIDERMIS MORPHOGENESIS GO; POSITIVE REGULATION OF SMOOTHENED SIGNALING PATHWAY GO.",
            "summary": "A likely causal atrial fibrillation GWAS variant, found inside a cardiac DHS as part of a much larger cardiac and musculoskeletal DHS domain"
        },
        {
            "query": "PLAG1 MOTIF enrichment @ 64kbp; Satellite REPEAT enrichment @ 1Mbp; HINFP1/3 MOTIF enrichment @ 256kbp; KLF/SP/2 MOTIF enrichment @ 16kbp; Ebox/CACGTG/1 MOTIF enrichment @ 4kbp; Mean corpuscular hemoglobin GWAS occurrence @ 1bp; NKD2 GENE; SLC12A7 GENE; AMMONIUM TRANSMEMBRANE TRANSPORT GO; MONOATOMIC ANION HOMEOSTASIS GO; POSITIVE REGULATION OF PROTEIN MATURATION GO",
            "summary": "A likely causal red blood cell GWAS variant, found inside a myeloid/erythroid DHS contained in an active enhancer element."
        },
        {
            "query": "Lymphoid DHS enrichment @ 16kbp; IRF/2 MOTIF enrichment @ 4kbp; NRF1 MOTIF enrichment @ 1Mbp; ZNF320 MOTIF enrichment @ 256kbp; SREBF1 MOTIF enrichment @ 64kbp; MECP2 motif occurrence @ 1bp; TFAP2/1 MOTIF occurrence @ 16bp; KLF/SP/2 MOTIF enrichment @ 1kbp; CCDC22 GENE; FOXP3 GENE; NEGATIVE REGULATION OF NF KAPPAB TRANSCRIPTION FACTOR ACTIVITY GO; NEGATIVE REGULATION OF DNA BINDING TRANSCRIPTION FACTOR ACTIVITY GO; REGULATION OF DNA BINDING TRANSCRIPTION FACTOR ACTIVITY GO",
            "summary": "Weak enhancer element harboring an AP-2 transcription factor motif, residing in a larger domain of interferon-regulatory factor (IRF) protein binding sites and lymphoid DHSs. Co-located with the FOXP3 gene, an important immune system regulator."
        }
    ]
    @outlines.prompt
    def summary_prompt(query, examples):
        """You are a genomics researcher who is an expert in the field of genomics, tasked with narrating genomic regions such that a short sentence captures most of the standout information.							
        You are given a query consisting of an term-wise description of a certain region in the human genome, that may include things like transcription factor motif hits (MOTIF), chromatin state calls (CS), DNaseI Hypersensitive Site component annotations (DHS), interspersed repeats and low complexity DNA sequences (REPEAT), and associated GWAS traits (GWAS).							
        All item types are observed at a certain genomic scale, ranging from single basepair (1bp) to Mega basepair (1Mbp), and are listed in the query in descending order of prominance, so make sure take that into account.							
        Additionally, you are provided information on any overlapping or nearby genes (GENE), and associated Gene Ontology genesets (GO), which may constitute important information.							
        Your task is to generate a helpful one-sentence summary of the query, essentially providing a useful narrative of the genomic region.							
        If any of the provided terms do not seem relevant according to literature or otherwise, feel free to skip them in the narrative.							

        Examples
        --------
        {% for example in examples %}
        Query: {{ example.query }}
        Summary: {{ example.summary }}
        {% endfor %}

        Task
        --------
        Query: {{ query }}
        Summary:
        """
    @outlines.prompt
    def rag_prompt(query, examples, abstracts):
        """You are a genomics researcher who is an expert in the field of genomics, tasked with narrating genomic regions such that a short sentence captures most of the standout information.							
        You are given a query consisting of an term-wise description of a certain region in the human genome, that may include things like transcription factor motif hits (MOTIF), chromatin state calls (CS), DNaseI Hypersensitive Site component annotations (DHS), interspersed repeats and low complexity DNA sequences (REPEAT), and associated GWAS traits (GWAS).							
        All item types are observed at a certain genomic scale, ranging from single basepair (1bp) to Mega basepair (1Mbp), and are listed in the query in descending order of prominance, so make sure take that into account.							
        Furthermore, items may refer to single instances ("occurrence"), or a local increase in the number of occurrences ("enrichment"). This is important information to use.
        Additionally, you are provided information on any overlapping or nearby genes (GENE), and associated Gene Ontology genesets (GO), which may constitute important information.							
        You also have access to titles and abstracts of research articles that may be relevant to the query, so make sure to use these for additional context and writing style.							
        Your task is to generate a helpful one-sentence summary of the query, essentially providing a useful narrative of the genomic region.							
        If any of the provided terms do not seem relevant according to literature or otherwise, feel free to skip them in the narrative.							

        Examples
        --------
        {% for example in examples %}
        Query: {{ example.query }}
        Summary: {{ example.summary }}
        {% endfor %}

        Abstracts
        --------

        {% for abstract in abstracts %}
        Title: {{ abstract.full_title }}
        Abstract: {{ abstract.abstract }}

        {% endfor %}
        Task
        --------

        Query: {{ query}}
        Summary:
        """
    return {
        "model": model,
        "table": table,
        "summary_prompt": summary_prompt,
        "rag_prompt": rag_prompt,
        "examples": examples,
        "default_instructions": default_instructions
    }

# initialize the model and its components on startup
model_components = start_engine()
print(f"🤖 model initialized", flush=True)


@pubmed_summary_bp.route("pubmed_summary", methods=["POST"])
def pubmed_summary():
    request_info = request.get_json()
    query = request_info["query"]
    # narration = request_info.get("narration", None)

    print("🔍🗣️ query", query, flush=True)

    start = time.time()
    print("MODEL DEBUG", model_components["model"], flush=True)

    try:
        embeddings = model_components["model"].encode(query, normalize_embeddings=True)
        print(f"Embeddings time: {time.time() - start}", flush=True)
    except Exception as e:
        print(f"Error during encoding: {str(e)}", flush=True)
        return json.dumps({"error": f"Error during encoding: {e}"})

    start = time.time()
    results = model_components["table"].search(embeddings).metric("cosine").limit(10).to_list()
    print(f"Results time: {time.time() - start}", flush=True)

    start = time.time()
    prompt =  model_components["rag_prompt"](query, model_components["examples"], results)
    print(f"Prompt time: {time.time() - start}", flush=True)

    start = time.time()
    model = models.openai("gpt-4o-mini", api_key=os.environ["OPENAI_API_KEY"])
    print(f"Model time: {time.time() - start}", flush=True)

    start = time.time()
    generator = generate.text(model)
    print(f"Generator time: {time.time() - start}", flush=True)
    
    # Remove vector field from results before returning
    for result in results:
        del result["vector"]

    start = time.time()
    summary = generator(prompt)
    print(f"Summary time: {time.time() - start}", flush=True)

    request_id = str(uuid.uuid4())

    # # Log a sample request
    # request_data = {
    #     "timestamp": datetime.datetime.now(),
    #     "request_id": request_id,
    #     "type": "rag_generate",
    #     "query": query,
    #     "summary": summary,
    #     "results": json.dumps(results),
    #     "narration": json.dumps(narration)
    # }
    # self.logger.log_request(request_data)
    # # Force flush to disk (normally happens automatically at buffer limit)
    # self.logger.flush()

    return json.dumps({"summary": summary, "results": results, "request_id": request_id})