"""Central configuration for pegrna-db."""
import os
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()

# Paths
BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
RAW_PAPERS_DIR = DATA_DIR / "raw_papers"
DATABASE_PATH = BASE_DIR / os.getenv("DATABASE_PATH", "data/pegrna_database.db")

# GitHub raw URL for the database (avoids HF Space LFS limits)
GITHUB_DB_URL = "https://media.githubusercontent.com/media/JoeyRiepsaame/pegrna-db/master/data/pegrna_database.db"


def _ensure_database():
    """Download database from GitHub if not present locally."""
    if DATABASE_PATH.exists():
        return
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    # Try gzip fallback first
    import gzip, shutil
    _db_gz = DATABASE_PATH.with_suffix(".db.gz")
    if _db_gz.exists():
        print(f"Decompressing {_db_gz} -> {DATABASE_PATH} ...")
        with gzip.open(_db_gz, "rb") as f_in, open(DATABASE_PATH, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        print(f"Done. Database size: {DATABASE_PATH.stat().st_size / 1024 / 1024:.1f} MB")
        return
    # Download from GitHub (LFS-hosted)
    try:
        import urllib.request
        print(f"Downloading database from GitHub ({GITHUB_DB_URL})...")
        urllib.request.urlretrieve(GITHUB_DB_URL, str(DATABASE_PATH))
        print(f"Downloaded: {DATABASE_PATH} ({DATABASE_PATH.stat().st_size / 1024 / 1024:.1f} MB)")
    except Exception as e:
        print(f"WARNING: Could not download database: {e}")


_ensure_database()

# Ensure directories exist
DATA_DIR.mkdir(exist_ok=True)
RAW_PAPERS_DIR.mkdir(exist_ok=True)

# NCBI / PubMed
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "")
NCBI_TOOL = os.getenv("NCBI_TOOL", "pegrna-db")

# Anthropic / Claude
ANTHROPIC_API_KEY = os.getenv("ANTHROPIC_API_KEY", "")

# OpenRouter (fallback LLM)
OPENROUTER_API_KEY = os.getenv("OPENROUTER_API_KEY", "")

# ChatGPT (fallback LLM)
CHATGPT_API_KEY = os.getenv("CHATGPT_API_KEY", "")

# PubMed search defaults
PUBMED_SEARCH_TERMS = [
    '"prime editing" AND (pegRNA OR epegRNA OR "prime editing guide RNA")',
]
PUBMED_MAX_RESULTS = 500

# PMC BioC API
PMC_BIOC_BASE = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json"
PMC_SUPP_BASE = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/FAIR-SMART/BioC_json"

# Unpaywall
UNPAYWALL_BASE = "https://api.unpaywall.org/v2"
UNPAYWALL_EMAIL = NCBI_EMAIL

# LLM extraction
# ClinVar
CLINVAR_FTP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
CLINVAR_DATA_DIR = DATA_DIR / "clinvar"

# Ensembl API
ENSEMBL_REST_SERVER = "https://rest.ensembl.org"
ENSEMBL_GENOMES_REST_SERVER = "https://rest.ensemblgenomes.org"
SPLICE_SITE_WINDOW = 2  # bp window around exon boundaries

LLM_MODEL = "claude-sonnet-4-6"
LLM_MAX_TOKENS = 8192

# Column name patterns for rule-based extraction (case-insensitive)
# Patterns are tried in order; put most specific first
SPACER_PATTERNS = [
    "spacer_sequence", "spacer sequence", "spacer", "guide_sequence", "guide",
    "protospacer", "protospacer_sequence", "protospacer_guide",
    "target_sequence", "target site sequence",
    "sgrna", "grna", "grna_sequence",
    "spacer_seq", "guide_seq", "pegrna_spacer",
    "prime_editing_pegrna_spacer_seq", "prime_editing_pegRNA_spacer_seq",
    "pegrna spacer", "pegRNA spacer",
    "n20", "20mer", "crrna",
]
PBS_PATTERNS = [
    "pbs_sequence", "pbs_seq", "primer_binding", "primer binding site",
    "pbs", "primer_binding_site", "primerbindingsite",
    "primer_binding_sequence", "pbs_dna", "pbs_oligo",
]
RTT_PATTERNS = [
    "rtt_sequence", "rtt_seq", "rt_template", "rt_seq",
    "reverse transcriptase template", "rt template", "rtt",
    "rt_template_sequence", "rtt_dna", "rttemplate",
    "reverse_transcriptase", "template_sequence", "rt_oligo",
    "template", "ha_rtt",
]
# Combined extension patterns (RTT+PBS in a single column)
EXTENSION_PATTERNS = [
    "extension_seq", "3' extension", "3'-extension", "extension sequence",
    "3' extension sequence", "3'-extension sequence",
    "pbs-rtt", "rtt+pbs", "rt+pbs", "pbs+rtt",
    "prime_editing_pegrna_extension_seq", "prime_editing_pegRNA_extension_seq",
    "pegrna 3' extension", "pegRNA 3' extension",
]
EFFICIENCY_PATTERNS = [
    "editing_efficiency", "efficiency", "edit_eff", "pe_efficiency",
    "percent_editing", "percentage_editing", "% editing", "editing (%)", "editing rate",
    "averageedited", "percentageedited",
]
EDIT_TYPE_PATTERNS = [
    "edit_type", "mutation_type", "edit type", "mutation type",
    "variant_type", "correction_type",
]
CELL_TYPE_PATTERNS = [
    "cell_type", "cell_line", "cell line", "cell type", "cells",
]
GENE_PATTERNS = [
    "target_gene", "target gene", "target_site", "target site",
    "gene", "target",
]
PE_VERSION_PATTERNS = [
    "pe_version", "prime_editor", "prime editor", "pe version",
]
