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
LLM_MODEL = "claude-sonnet-4-5-20250929"
LLM_MAX_TOKENS = 4096

# Column name patterns for rule-based extraction (case-insensitive)
SPACER_PATTERNS = [
    "spacer", "guide", "protospacer", "target_sequence", "guide_sequence",
    "sgrna", "grna", "spacer_sequence",
]
PBS_PATTERNS = [
    "pbs", "primer_binding", "primer binding site", "pbs_sequence",
    "pbs_length", "pbs_len",
]
RTT_PATTERNS = [
    "rtt", "rt_template", "reverse transcriptase template", "rtt_sequence",
    "rtt_length", "rtt_len", "rt template", "template",
]
EFFICIENCY_PATTERNS = [
    "efficiency", "editing_efficiency", "edit_eff", "pe_efficiency",
    "percent_editing", "% editing", "editing (%)", "editing rate",
]
EDIT_TYPE_PATTERNS = [
    "edit_type", "mutation_type", "edit type", "mutation type",
    "variant_type", "substitution", "insertion", "deletion",
]
CELL_TYPE_PATTERNS = [
    "cell_type", "cell_line", "cell line", "cell type", "cells",
]
GENE_PATTERNS = [
    "gene", "target_gene", "locus", "target_locus", "target gene",
    "target site", "target_site",
]
PE_VERSION_PATTERNS = [
    "pe_version", "prime_editor", "editor", "pe version", "pe2", "pe3",
    "pe4", "pe5", "pe7", "prime editor",
]
