from flask import Blueprint, request, g
import json
import os
import datetime
import time

os.environ["TOKENIZERS_PARALLELISM"] = "false"
os.environ["HF_HOME"]="~/.cache/huggingface/datasets"

from sentence_transformers import SentenceTransformer
import outlines
from outlines import models, generate
import torch
import uuid
import pyarrow as pa
from logger import ParquetRequestLogger
from constants import MODEL_ID, LANCEDB_URI, LOGGING_DIR

api_key = os.environ["OPENAI_API_KEY"]

# Create a Blueprint
pubmed_summary_bp = Blueprint('pubmed_summary_bp', __name__)


_db_connection = None
def get_db_connection():
    """
    Collects lancedb connection for current worker. If the connection is not already open, 
    opens a new connection.
    """

    import lancedb  # import here to ensure that the module is available for the current worker

    # for current worker, checks to see if the DB connection is already open, and if not, opens a new connection
    global _db_connection
    if _db_connection is None:
        print("Initializing new DB connection", flush=True)
        _db_connection = lancedb.connect(LANCEDB_URI)
    return _db_connection


_model_components = None
def get_model_components():
    """
    Collects model components for the current worker. If the model is not already loaded, 
    loads a new model instance for worker.
    """

    # for current worker, checks to see if the model is already loaded, and if not, loads a new model instance
    global _model_components
    if _model_components is None:
        device = torch.device("cpu")
        print("🥶 cold starting inference", flush=True)
        start = time.monotonic_ns()

        # loads model
        model = SentenceTransformer(MODEL_ID, trust_remote_code=True, device=device)
        model.eval()
        duration_s = (time.monotonic_ns() - start) / 1e9
        print(f"🏎️ engine started in {duration_s:.0f}s", flush=True)

        # establishes lancedb connection and loads table
        db = get_db_connection() 
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

        # # Create logger with custom schema (optional)
        # request_schema = pa.schema([
        #     ('timestamp', pa.timestamp('us')),
        #     ('request_id', pa.string()),
        #     ('type', pa.string()),
        #     ('query', pa.string()),
        #     ('prompt', pa.string()),
        #     ('summary', pa.string()),
        #     ('results', pa.string()),
        #     ('narration', pa.string())
        # ])
        # request_logger = ParquetRequestLogger(log_dir=os.path.join(LOGGING_DIR, "requests"), schema=request_schema)

        _model_components = {
            "model": model,
            "table": table,
            "summary_prompt": summary_prompt,
            "rag_prompt": rag_prompt,
            "examples": examples,
            "default_instructions": default_instructions,
            # "request_logger": request_logger
        }
    return _model_components


def initialize_request_logger():
    # Create logger with custom schema (optional) shared by all workers
    request_schema = pa.schema([
        ('timestamp', pa.timestamp('us')),
        ('request_id', pa.string()),
        ('type', pa.string()),
        ('query', pa.string()),
        ('prompt', pa.string()),
        ('summary', pa.string()),
        ('results', pa.string()),
        ('narration', pa.string())
    ])
    request_logger = ParquetRequestLogger(log_dir=os.path.join(LOGGING_DIR, "requests"), schema=request_schema)
    return request_logger
# logger shared by all workers
request_logger = initialize_request_logger()



def get_summary(query, st_model, table, rag_prompt, examples, api_key):
    """
    Generates a text summary for a given query using the provided sentence transformers model, 
    lancedb table, rag_prompt, and examples. Returns sentence summary and search results.
    """

    start = time.time()
    # collect embeddings for the query
    embeddings = st_model.encode(query, normalize_embeddings=True)
    print(f"\tEmbeddings time: {time.time() - start}", flush=True)
    
    start = time.time()
    # search table
    results = table.search(embeddings).metric("cosine").limit(10).to_list()
    print(f"\tResults time: {time.time() - start}", flush=True)

    start = time.time()
    # generate prompt
    prompt =  rag_prompt(query, examples, results)
    print(f"\tPrompt time: {time.time() - start}", flush=True)

    start = time.time()
    # load openAI model
    model = models.openai("gpt-4o-mini", api_key=api_key)
    print(f"\tModel time: {time.time() - start}", flush=True)

    start = time.time()
    # create generator
    generator = generate.text(model)
    print(f"\tGenerator time: {time.time() - start}", flush=True)

    # Remove vector field from results before returning
    for result in results:
        del result["vector"]

    start = time.time()
    # feed prompt to generator
    summary = generator(prompt)
    print(f"\tSummary time: {time.time() - start}", flush=True)

    return summary, results, prompt


def create_custom_prompt(prompt_template):
        """
        Defines a custom prompt function for the summary generation task using a user-provided prompt template.
        """
        def custom_prompt(query, examples, abstracts):
            return ""
        custom_prompt.__doc__ = prompt_template
        return outlines.prompt(custom_prompt)


@pubmed_summary_bp.route("pubmed_summary", methods=["POST"])
def pubmed_summary():
    """
    Endpoint for generating a summary for a given query
    """

    request_info = request.get_json()
    query = request_info["query"]  # 'Cancer / epithelial DHS enrichment @ 64kbp'
    narration = request_info.get("narration", None)
    prompt = request_info.get("prompt", None)
    print("🔍🗣️ query", query, flush=True)
    
    # collects model components for the current worker
    model_components = get_model_components()

    rag_prompt = model_components["rag_prompt"] if prompt is None else create_custom_prompt(prompt)

    # generates summary and results
    summary, results, used_prompt = get_summary(query, model_components["model"], model_components["table"], rag_prompt, model_components["examples"], api_key)
    request_id = str(uuid.uuid4())

    # log request
    request_data = {
        "timestamp": datetime.datetime.now(),
        "request_id": request_id,
        "type": "rag_generate",
        "query": query,
        "prompt": used_prompt,
        "summary": summary,
        "results": json.dumps(results),
        "narration": json.dumps(narration)
    }
    request_logger.log_request(request_data)
    # Force flush to disk (normally happens automatically at buffer limit)
    request_logger.flush()

    return json.dumps({"summary": summary, "results": results, "request_id": request_id})