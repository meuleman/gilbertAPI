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
                "query": "EWSR1/FLI1 MOTIF domain @ 16kbp; Stromal B DHS domain @ 1Mbp; Atrial fibrillation GWAS occurrence @ 1bp; Musculoskeletal DHS domain @ 64kbp; Quiescent/Low CS occurrence @ 256bp; Cardiac DHS domain @ 256kbp; Quiescent/Low CS domain @ 1kbp; Cardiac DHS occurrence @ 64bp; NTMT2 GENE_ADJ; GORAB GENE_ADJ; N TERMINAL PROTEIN AMINO ACID MODIFICATION GO; EPIDERMIS MORPHOGENESIS GO; POSITIVE REGULATION OF SMOOTHENED SIGNALING PATHWAY GO.",
                "summary": "A likely causal atrial fibrillation GWAS variant, found inside a cardiac DHS as part of much larger cardiac and musculoskeletal DHS domains"
            },
            {
                "query": "PLAG1 MOTIF domain @ 64kbp; Satellite REPEAT domain @ 1Mbp; HINFP1/3 MOTIF domain @ 256kbp; KLF/SP/2 MOTIF domain @ 16kbp; Ebox/CACGTG/1 MOTIF domain @ 4kbp; Mean corpuscular hemoglobin GWAS occurrence @ 1bp; Myeloid / erythroid DHS occurrence @ 256bp; NKD2 GENE_ADJ; SLC12A7 GENE_OVL; AMMONIUM TRANSMEMBRANE TRANSPORT GO; MONOATOMIC ANION HOMEOSTASIS GO; POSITIVE REGULATION OF PROTEIN MATURATION GO",
                "summary": "A likely causal red blood cell GWAS variant, found inside a myeloid/erythroid DHS."
            },
            {
                "query": "Lymphoid DHS domain @ 16kbp; IRF/2 MOTIF domain @ 4kbp; NRF1 MOTIF domain @ 1Mbp; ZNF320 MOTIF domain @ 256kbp; SREBF1 MOTIF domain @ 64kbp; KLF/SP/2 MOTIF domain @ 1kbp; KLF/SP/2 MOTIF occurrence @ 16bp; CCDC22 GENE_ADJ; FOXP3 GENE_OVL; NEGATIVE REGULATION OF NF KAPPAB TRANSCRIPTION FACTOR ACTIVITY GO; NEGATIVE REGULATION OF DNA BINDING TRANSCRIPTION FACTOR ACTIVITY GO; REGULATION OF DNA BINDING TRANSCRIPTION FACTOR ACTIVITY GO",
                "summary": "Region residing in a larger domain of interferon-regulatory factor (IRF) protein binding sites and lymphoid DHSs, overlapping the FOXP3 gene, an important immune system regulator."
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
            """You are an expert genomics researcher, tasked with narrating genomic regions such that a single short sentence captures the most important information.
            You are given a query consisting of a term-wise description of a certain region of interest in the human genome.

            These terms may include things like chromatin state calls (CS), DNase I Hypersensitive Site annotations (DHS), transcription factor motif hits (MOTIF), interspersed repeats and low complexity DNA sequences (REPEAT), and Genome-Wide Association Study traits (GWAS).
            All such terms are observed at a certain genomic scale, ranging from a single basepair (1bp) to a million basepair (1Mbp).
            They are listed in the query in descending order of prominence, so make sure to take that into account in prioritizing the information to use in your narration.
            Furthermore, the genomic region of interest may directly overlap an observed term ('occurrence'), or may overlap a larger region with an abundance of that term ('domain') in which there is not necessarily a direct overlap with a single instance of the term. This is an important distinction.

            Additionally, you are provided information on any genes that directly overlap (GENE_OVL) or are adjacent to (GENE_ADJ) the region.
            To aid in functional narration, you are also provided with Gene Ontology genesets (GO) associated with the region, which may constitute important information in combination with all of the above.

            You also have access to titles and abstracts of research articles that may be relevant to the query, so make sure to use these for additional context and writing style.							

            Your task is to generate a helpful one-sentence summary of the query, providing a useful narrative of the genomic region.
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
        region_set_examples = [
            # {
            #     "query": "Cardiac DHS domain @ 16kbp: 81%; Cardiac DHS domain @ 64kbp: 19%; HD/18 MOTIF domain @ 1Mbp: 11%; HD/18 MOTIF domain @ 256kbp: 10%",
            #     "summary": "A set of cardiac DHS domains most commonly at the 16kbp scale, featuring significant motif enrichment for HD/18 transcription factors."
            # }
            {
                "query": "Eosinophil count GWAS  @ 1bp 74%; Lymphoid DHS domain @ 64kbp 30%; Myeloid / erythroid DHS domain @ 64kbp 24%; KLF/SP/2 MOTIF domain @ 16kbp 23%; MECP2 MOTIF domain @ 256kbp 20%; ZNF384/2 MOTIF domain @ 1Mbp 18%; Weak transcription CS domain @ 4kbp 15%; Weak Repressed PolyComb CS domain @ 1kbp 13%; Quiescent/Low CS domain @ 4kbp 13%; Placental / trophoblast DHS domain @ 1Mbp 12%; GO GOBP CELL ACTIVATION; GO GOBP MONONUCLEAR CELL DIFFERENTIATION; GO GOBP LEUKOCYTE DIFFERENTIATION; GENE_OVL IL18R1; GENE_ADJ IL18RAP; GENE_ADJ IL5; GENE_ADJ IL13; GENE_OVL CDK6; GENE_ADJ TSC1; GENE_OVL CEBPE; GENE_OVL BCL2; GENE_ADJ CSF2; GENE_ADJ CCR7; GENE_ADJ SMARCE1; GENE_ADJ TSLP; GENE_OVL SMAD3; GENE_OVL TNFRSF1B; GENE_OVL GATA3; GENE_OVL ICOSLG; GENE_OVL IL17RA; GENE_OVL IRF4; GENE_ADJ BATF; GENE_ADJ LRRC32; GENE_ADJ PLCG2; GENE_OVL PTPRC; GENE_ADJ GATA2; GENE_OVL BCL3; GENE_ADJ IL33; GENE_ADJ FOXP1; GENE_OVL JAK1; GENE_OVL HMGB1; GENE_ADJ GPR183; GENE_ADJ ID2; GENE_OVL PRG3; GENE_OVL JAK2; GENE_ADJ ITGB8; GENE_OVL CLC; GENE_OVL IKZF1",
                "summary": "This set contains regions with single nucleotide polymorphisms significantly associated with eosinophil count identified through GWAS, DHS domains indicative of regulatory elements involved in lymphoid and myeloid cell activation and differentiation, enrichment for multiple transcription factor motifs, and is linked with pathways related to immune response through overlapping genes such as IL5, IL13, and CEBPE."
            },
            {
                "query": "Cardiac DHS domain @ 16kbp 100%; HD/18 MOTIF domain @ 1Mbp 28%; EWSR1/FLI1 MOTIF domain @ 64kbp 18%; KLF/SP/2 MOTIF domain @ 64kbp 11%; Placental / trophoblast DHS domain @ 1Mbp 8%; Stromal B DHS domain @ 1Mbp 8%; ZFN121 MOTIF domain @ 256kbp 7%; HD/2 MOTIF domain @ 1Mbp 5%; Vascular / endothelial DHS domain @ 1Mbp 5%; Digestive DHS domain @ 1Mbp 5%; GO GOBP STRIATED MUSCLE TISSUE DEVELOPMENT; GO GOBP HEART DEVELOPMENT; GO GOBP CARDIAC CHAMBER DEVELOPMENT; GENE_OVL SLC8A1; GENE_ADJ PDLIM3; GENE_OVL SORBS2; GENE_OVL PLN; GENE_OVL BMP5; GENE_OVL AKAP6; GENE_OVL RBM24; GENE_OVL TBX20; GENE_OVL POPDC2; GENE_OVL GATA6; GENE_OVL LDB3; GENE_ADJ BMPR1A; GENE_ADJ GATA5; GENE_OVL ACTN2; GENE_OVL GATA4; GENE_OVL ZFPM2; GENE_ADJ GJA5; GENE_OVL TNNT2; GENE_OVL HAND2; GENE_OVL CACNA1C; GENE_OVL MIB1; GENE_OVL HEY2; GENE_OVL MYH6; GENE_OVL MYH7; GENE_OVL NDST1; GENE_OVL NKX2-5; GENE_ADJ TBX5; GENE_ADJ TBX3; GENE_OVL MYOCD; GENE_ADJ MTPN; GENE_OVL PROX1; GENE_ADJ SMYD2; GENE_ADJ FZD1; GENE_ADJ TMEM65; GENE_ADJ ANKRD1; GENE_OVL NEBL",
                "summary": "The genomic regions predominantly feature cardiac DHS domains at 16kbp, enrichments for various transcription factor motifs, and significant overlaps with crucial cardiac genes like NKX2-5 and MYH6, reflecting their role in striated muscle tissue and cardiac chamber development."
            },
            {
                "query": "Active Enhancer 1 CS domain @ 256bp 100%; Active Enhancer 1 CS occurrence @ 64bp 98%; CTCF MOTIF occurrence @ 16bp 89%; KLF/SP/2 MOTIF domain @ 16kbp 18%; Primitive / embryonic DHS occurrence @ 64bp 17%; ZNF384/2 MOTIF domain @ 1Mbp 14%; MECP2 MOTIF domain @ 256kbp 13%; ZFN121 MOTIF domain @ 1Mbp 9%; ZNF768 MOTIF domain @ 1Mbp 9%; Placental / trophoblast DHS domain @ 64kbp 9%; GO GOBP POSITIVE REGULATION OF NUCLEAR TRANSCRIBED MRNA CATABOLIC PROCESS DEADENYLATION DEPENDENT DECAY; GO GOBP REGULATION OF NUCLEAR TRANSCRIBED MRNA CATABOLIC PROCESS DEADENYLATION DEPENDENT DECAY; GO GOBP ARTERY DEVELOPMENT; GENE_ADJ TNRC6C; GENE_ADJ TOB1; GENE_ADJ ZFP36L2; GENE_ADJ NAGLU; GENE_OVL SMAD7; GENE_OVL ZMIZ1; GENE_ADJ ZFP36L1; GENE_ADJ PLXND1; GENE_OVL SMAD6; GENE_ADJ HES1; GENE_OVL SUFU",
                "summary": "These regions comprise highly active enhancer chromatin state domains and significant CTCF motif occurrences, with adjacent genes involved in mRNA regulation and developmental processes, suggesting a complex interplay of transcriptional control and enhancer dynamics."
            }
        ]
        @outlines.prompt
        def region_set_prompt(query, examples, abstracts):
            """You are an expert genomics researcher, tasked with narrating a set of genomic regions such that a single short sentence captures the most important information about them.
            You are given a query consisting of a term-wise description of a set of regions across the human genome.

            These terms may include things like chromatin state calls (CS), DNase I Hypersensitive Site annotations (DHS), transcription factor motif hits (MOTIF), interspersed repeats and low complexity DNA sequences (REPEAT), and Genome-Wide Association Study traits (GWAS).
            All such terms are observed at a certain genomic scale, ranging from a single basepair (1bp) to a million basepair (1Mbp).
            They are listed in the query in descending order of prominence, with the percent of regions in the set containing each term also listed, so make sure to take that into account in prioritizing the information to use in your narration.
            Furthermore, the genomic regions of interest may directly overlap an observed term ('occurrence'), or may overlap a larger region with an abundance of that term ('domain') in which there is not necessarily a direct overlap with a single instance of the term. This is an important distinction.

            Additionally, you are provided information on any genes that directly overlap (GENE_OVL) or are adjacent to (GENE_ADJ) the regions in the set.
            To aid in the functional narration, you are also provided with Gene Ontology genesets (GO) associated with the regions, including corresponding p-values that quantify the strength of the association, which may constitute important information in combination with all of the above.

            You also have access to titles and abstracts of research articles that may be relevant to the query, so make sure to use these for additional context and writing style.							

            Your task is to generate a helpful one-sentence summary of the query, providing a useful narrative of the genomic regions.
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
            "region_set_prompt": region_set_prompt,
            "region_set_examples": region_set_examples,
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


def initialize_feedback_logger():
    # Create logger with custom schema (optional) shared by all workers
    feedback_schema = pa.schema([
        ('timestamp', pa.timestamp('us')),
        ('request_id', pa.string()),
        ('feedback', pa.string())
    ])
    feedback_logger = ParquetRequestLogger(log_dir=os.path.join(LOGGING_DIR, "feedback"), schema=feedback_schema)
    return feedback_logger
# logger shared by all workers
feedback_logger = initialize_feedback_logger()



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
    prompt = rag_prompt(query, examples, results)
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


@pubmed_summary_bp.route("pubmed_region_set_summary", methods=["POST"])
def pubmed_region_set_summary():
    """
    Endpoint for generating a summary for a set of regions
    """

    request_info = request.get_json()
    query = request_info["query"]  # 'Cancer / epithelial DHS enrichment @ 64kbp | ...'
    # narration = request_info.get("narration", None)
    print("🔍🗣️ query", query, flush=True)
    
    # collects model components for the current worker
    model_components = get_model_components()

    region_set_prompt = model_components["region_set_prompt"]

    # generates summary and results
    summary, results, used_prompt = get_summary(query, model_components["model"], model_components["table"], region_set_prompt, model_components["region_set_examples"], api_key)
    request_id = str(uuid.uuid4())

    # # log request
    # request_data = {
    #     "timestamp": datetime.datetime.now(),
    #     "request_id": request_id,
    #     "type": "rag_generate",
    #     "query": query,
    #     "prompt": used_prompt,
    #     "summary": summary,
    #     "results": json.dumps(results),
    #     "narration": json.dumps(narration)
    # }
    # request_logger.log_request(request_data)
    # # Force flush to disk (normally happens automatically at buffer limit)
    # request_logger.flush()

    return json.dumps({"summary": summary, "results": results, "request_id": request_id})


@pubmed_summary_bp.route("feedback", methods=["POST"])
def feedback():
    """
    Endpoint for providing feedback on the model
    """

    request_info = request.get_json()

    feedback_logger.log_request({
        "timestamp": datetime.datetime.now(),
        "request_id": request_info["request_id"],
        "feedback": request_info["feedback"],
    })
    # Force flush to disk (normally happens automatically at buffer limit)
    feedback_logger.flush()

    return json.dumps({"status": "success"})


@pubmed_summary_bp.route("get_single_region_prompt", methods=["POST"])
def get_single_region_prompt():
    """
    Endpoint for generating a summary for a set of regions
    """

    # collects model components for the current worker
    model_components = get_model_components()

    return json.dumps({"prompt": getattr(model_components["rag_prompt"].__doc__, "template", None)})



@pubmed_summary_bp.route("get_region_set_prompt", methods=["POST"])
def get_region_set_prompt():
    """
    Endpoint for generating a summary for a set of regions
    """

    # collects model components for the current worker
    model_components = get_model_components()

    return json.dumps({"prompt": getattr(model_components["region_set_prompt"].__doc__, "template", None)})