"""Ensembl REST API client for gene structure retrieval and caching."""
import json
import time
from datetime import datetime
from typing import Optional

import requests
from rich.console import Console
from sqlalchemy.orm import Session

from database.models import GeneStructure

console = Console()

# Map target_organism values -> (ensembl_species, api_server)
ORGANISM_MAP = {
    "Homo sapiens":            ("homo_sapiens",            "https://rest.ensembl.org"),
    "Human":                   ("homo_sapiens",            "https://rest.ensembl.org"),
    "human":                   ("homo_sapiens",            "https://rest.ensembl.org"),
    "Mus musculus":            ("mus_musculus",             "https://rest.ensembl.org"),
    "Mouse":                   ("mus_musculus",             "https://rest.ensembl.org"),
    "mouse":                   ("mus_musculus",             "https://rest.ensembl.org"),
    "Danio rerio":             ("danio_rerio",             "https://rest.ensembl.org"),
    "Zebrafish":               ("danio_rerio",             "https://rest.ensembl.org"),
    "zebrafish":               ("danio_rerio",             "https://rest.ensembl.org"),
    "Rattus norvegicus":       ("rattus_norvegicus",       "https://rest.ensembl.org"),
    "Sus scrofa":              ("sus_scrofa",              "https://rest.ensembl.org"),
    "Drosophila melanogaster": ("drosophila_melanogaster", "https://rest.ensembl.org"),
    "Oryza sativa":            ("oryza_sativa",            "https://rest.ensemblgenomes.org"),
    "rice":                    ("oryza_sativa",            "https://rest.ensemblgenomes.org"),
    "Rice":                    ("oryza_sativa",            "https://rest.ensemblgenomes.org"),
    "Arabidopsis thaliana":    ("arabidopsis_thaliana",    "https://rest.ensemblgenomes.org"),
    "Glycine max":             ("glycine_max",             "https://rest.ensemblgenomes.org"),
    "Zea mays":                ("zea_mays",                "https://rest.ensemblgenomes.org"),
    "Triticum aestivum":       ("triticum_aestivum",       "https://rest.ensemblgenomes.org"),
    "Nicotiana benthamiana":   ("nicotiana_benthamiana",   "https://rest.ensemblgenomes.org"),
    "Solanum lycopersicum":    ("solanum_lycopersicum",    "https://rest.ensemblgenomes.org"),
    "Physcomitrium patens":    ("physcomitrium_patens",    "https://rest.ensemblgenomes.org"),
}

RATE_LIMIT_DELAY = 0.07  # ~14 req/s
MAX_RETRIES = 3


def _request_with_retry(url: str, method: str = "GET", json_data=None,
                         headers=None, retries: int = MAX_RETRIES) -> Optional[requests.Response]:
    """Make an HTTP request with exponential backoff on rate limit."""
    if headers is None:
        headers = {"Content-Type": "application/json", "Accept": "application/json"}

    for attempt in range(retries):
        time.sleep(RATE_LIMIT_DELAY)
        try:
            if method == "POST":
                resp = requests.post(url, json=json_data, headers=headers, timeout=30)
            else:
                resp = requests.get(url, headers=headers, timeout=30)

            if resp.status_code == 200:
                return resp
            elif resp.status_code == 429:
                wait = 2 ** (attempt + 1)
                console.print(f"[yellow]Rate limited, waiting {wait}s...[/yellow]")
                time.sleep(wait)
            elif resp.status_code == 400:
                return None  # Bad request (gene not found etc.)
            else:
                console.print(f"[yellow]Ensembl API {resp.status_code}: {resp.text[:100]}[/yellow]")
                return None
        except requests.RequestException as e:
            console.print(f"[yellow]Request error: {e}[/yellow]")
            if attempt < retries - 1:
                time.sleep(2 ** (attempt + 1))

    return None


def lookup_genes_batch(species: str, symbols: list[str],
                       server: str) -> dict:
    """Batch lookup gene info from Ensembl. Returns {symbol: gene_data}."""
    url = f"{server}/lookup/symbol/{species}"
    result = {}

    # Ensembl accepts up to 1000 symbols per POST
    for i in range(0, len(symbols), 500):
        batch = symbols[i:i + 500]
        resp = _request_with_retry(
            url,
            method="POST",
            json_data={"symbols": batch, "expand": 1},
        )
        if resp:
            data = resp.json()
            for sym, gene_data in data.items():
                if isinstance(gene_data, dict) and gene_data.get("id"):
                    result[sym] = gene_data
        else:
            # Try individual lookups for the batch as fallback
            for sym in batch:
                single_url = f"{server}/lookup/symbol/{species}/{sym}?expand=1"
                r = _request_with_retry(single_url)
                if r:
                    result[sym] = r.json()

    return result


def get_canonical_transcript(gene_data: dict) -> Optional[dict]:
    """Extract the canonical transcript from gene lookup response."""
    transcripts = gene_data.get("Transcript", [])
    if not transcripts:
        return None

    # Find canonical
    for t in transcripts:
        if t.get("is_canonical") == 1:
            return t

    # Fallback: longest protein-coding transcript
    coding = [t for t in transcripts if t.get("biotype") == "protein_coding"]
    if coding:
        return max(coding, key=lambda t: abs(t.get("end", 0) - t.get("start", 0)))

    # Last resort: longest transcript
    return max(transcripts, key=lambda t: abs(t.get("end", 0) - t.get("start", 0)))


def get_exon_ranges(transcript: dict) -> list[tuple[int, int]]:
    """Extract sorted exon (start, end) ranges from transcript data."""
    exons = transcript.get("Exon", [])
    if not exons:
        return []

    ranges = [(e["start"], e["end"]) for e in exons if "start" in e and "end" in e]
    ranges.sort(key=lambda x: x[0])
    return ranges


def fetch_gene_sequence(ensembl_gene_id: str, server: str) -> Optional[str]:
    """Fetch the full genomic DNA sequence of a gene."""
    url = f"{server}/sequence/id/{ensembl_gene_id}?type=genomic"
    resp = _request_with_retry(url)
    if resp:
        data = resp.json()
        return data.get("seq")
    return None


def get_or_fetch_gene_structure(
    session: Session, gene_symbol: str, organism: str
) -> Optional[GeneStructure]:
    """Get gene structure from cache or fetch from Ensembl API."""
    # Check cache
    cached = (
        session.query(GeneStructure)
        .filter_by(gene_symbol=gene_symbol, organism=organism)
        .first()
    )
    if cached and cached.fetch_status in ("success", "not_found"):
        return cached if cached.fetch_status == "success" else None

    # Resolve organism
    org_info = ORGANISM_MAP.get(organism)
    if not org_info:
        # Try case-insensitive match
        for k, v in ORGANISM_MAP.items():
            if k.lower() == organism.lower():
                org_info = v
                break
    if not org_info:
        return None

    species, server = org_info

    # Lookup gene
    genes = lookup_genes_batch(species, [gene_symbol], server)
    gene_data = genes.get(gene_symbol)

    if not gene_data:
        # Store not_found to avoid re-fetching
        gs = GeneStructure(
            gene_symbol=gene_symbol,
            organism=organism,
            fetch_status="not_found",
            fetch_date=datetime.utcnow(),
        )
        session.add(gs)
        session.flush()
        return None

    transcript = get_canonical_transcript(gene_data)
    if not transcript:
        gs = GeneStructure(
            gene_symbol=gene_symbol,
            organism=organism,
            ensembl_gene_id=gene_data.get("id"),
            fetch_status="not_found",
            fetch_date=datetime.utcnow(),
        )
        session.add(gs)
        session.flush()
        return None

    exon_ranges = get_exon_ranges(transcript)

    # Fetch genomic sequence
    gene_seq = fetch_gene_sequence(gene_data["id"], server)

    gs = GeneStructure(
        gene_symbol=gene_symbol,
        organism=organism,
        ensembl_gene_id=gene_data.get("id"),
        ensembl_transcript_id=transcript.get("id"),
        chromosome=str(gene_data.get("seq_region_name", "")),
        strand=gene_data.get("strand"),
        gene_start=gene_data.get("start"),
        gene_end=gene_data.get("end"),
        exon_coordinates=json.dumps(exon_ranges),
        gene_sequence=gene_seq,
        fetch_status="success" if gene_seq and exon_ranges else "error",
        fetch_date=datetime.utcnow(),
    )
    session.add(gs)
    session.flush()

    return gs if gs.fetch_status == "success" else None


def batch_fetch_gene_structures(
    session: Session, gene_symbols: list[str], organism: str
) -> dict[str, GeneStructure]:
    """Fetch gene structures for multiple genes efficiently. Returns {symbol: GeneStructure}."""
    result = {}

    # Check cache first
    cached = (
        session.query(GeneStructure)
        .filter(
            GeneStructure.gene_symbol.in_(gene_symbols),
            GeneStructure.organism == organism,
        )
        .all()
    )
    cached_map = {gs.gene_symbol: gs for gs in cached}

    for sym in gene_symbols:
        if sym in cached_map:
            gs = cached_map[sym]
            if gs.fetch_status == "success":
                result[sym] = gs
            continue

    # Fetch missing genes
    missing = [s for s in gene_symbols if s not in cached_map]
    if not missing:
        return result

    org_info = ORGANISM_MAP.get(organism)
    if not org_info:
        for k, v in ORGANISM_MAP.items():
            if k.lower() == organism.lower():
                org_info = v
                break
    if not org_info:
        return result

    species, server = org_info

    # Batch lookup
    genes_data = lookup_genes_batch(species, missing, server)

    for sym in missing:
        gene_data = genes_data.get(sym)
        if not gene_data:
            gs = GeneStructure(
                gene_symbol=sym, organism=organism,
                fetch_status="not_found", fetch_date=datetime.utcnow(),
            )
            session.add(gs)
            continue

        transcript = get_canonical_transcript(gene_data)
        if not transcript:
            gs = GeneStructure(
                gene_symbol=sym, organism=organism,
                ensembl_gene_id=gene_data.get("id"),
                fetch_status="not_found", fetch_date=datetime.utcnow(),
            )
            session.add(gs)
            continue

        exon_ranges = get_exon_ranges(transcript)
        gene_seq = fetch_gene_sequence(gene_data["id"], server)

        gs = GeneStructure(
            gene_symbol=sym, organism=organism,
            ensembl_gene_id=gene_data.get("id"),
            ensembl_transcript_id=transcript.get("id"),
            chromosome=str(gene_data.get("seq_region_name", "")),
            strand=gene_data.get("strand"),
            gene_start=gene_data.get("start"),
            gene_end=gene_data.get("end"),
            exon_coordinates=json.dumps(exon_ranges),
            gene_sequence=gene_seq,
            fetch_status="success" if gene_seq and exon_ranges else "error",
            fetch_date=datetime.utcnow(),
        )
        session.add(gs)
        if gs.fetch_status == "success":
            result[sym] = gs

    session.flush()
    return result
