from flask import Flask, request, abort, Response
from flask_cors import CORS
import sys
from git import Repo


# CSN
from csns import csns_bp
from data_filtering import data_filtering_bp
from gwas import gwas_bp
from geneset_enrichment import geneset_enrichment_bp
from filtering_without_order import filtering_without_order_bp
from region_set_enrichment import region_set_enrichment_bp
from genes import genes_bp
from pubmed_summary import pubmed_summary_bp
from umap import umap_bp

app = Flask(__name__)
CORS(app)
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024 # 16 MB


current_name = (__file__)
prod_name = '/data/narrativeGenomics/NarrativeGenomics/region_narration/region_narration_app/region-narration-flask/region_narration.py'

if current_name == prod_name:
    server_root_dir = '/data/narrativeGenomics/NarrativeGenomics/region_narration'

    RN_VERSION = '0.0.1'
    RN_VERSION_STR = "{}/{}".format(Repo(server_root_dir, search_parent_directories=True).head.ref.name, RN_VERSION)

else:
    RN_VERSION_STR = ''


@app.route('/', methods=['GET'])
def about():
    python_version = sys.version
    return "<html><body><div>region-narration ({})</div><div>Python {}</div></body></html>".format(RN_VERSION_STR, python_version)


#############CSN############
app.register_blueprint(csns_bp, url_prefix='/api/csns') # Register the Blueprint at /api/csns
##########Filtering#########
app.register_blueprint(data_filtering_bp, url_prefix='/api/dataFiltering') # Register the Blueprint at /api/dataFiltering
############GWAS############
app.register_blueprint(gwas_bp, url_prefix='/api/gwas') # Register the Blueprint at /api/gwas
#####Geneset Enrichment#####
app.register_blueprint(geneset_enrichment_bp, url_prefix='/api/genesetEnrichment') # Register the Blueprint at /api/genesetEnrichment
##Filtering Without Order###
app.register_blueprint(filtering_without_order_bp, url_prefix='/api/filteringWithoutOrder') # Register the Blueprint at /api/filteringWithoutOrder
###Region Set Enrichment####
app.register_blueprint(region_set_enrichment_bp, url_prefix='/api/regionSetEnrichment') # Register the Blueprint at /api/regionSetEnrichment
###########Genes############
app.register_blueprint(genes_bp, url_prefix='/api/genes') # Register the Blueprint at /api/genes
####Pubmed Summary Model####
app.register_blueprint(pubmed_summary_bp, url_prefix='/api/pubmedSummary') # Register the Blueprint at /api/pubmedSummary
############UMAP############
app.register_blueprint(umap_bp, url_prefix='/api/umap') # Register the Blueprint at /api/umap

if __name__ == '__main__':
    app.run(host='0.0.0.0')
