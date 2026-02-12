"""Manual PMID/DOI input handling."""
import re
from pathlib import Path
from rich.console import Console

console = Console()

# Patterns
PMID_PATTERN = re.compile(r"^\d{6,9}$")
PMCID_PATTERN = re.compile(r"^PMC\d+$", re.IGNORECASE)
DOI_PATTERN = re.compile(r"^10\.\d{4,}/\S+$")


def parse_identifier(identifier: str) -> dict:
    """Parse a paper identifier (PMID, PMCID, DOI, or URL).

    Returns dict with keys: type, value (normalized).
    """
    identifier = identifier.strip()

    # Extract PMID from PubMed URL
    pmid_url = re.search(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)", identifier)
    if pmid_url:
        return {"type": "pmid", "value": pmid_url.group(1)}

    # Extract PMCID from PMC URL
    pmc_url = re.search(r"pmc\.ncbi\.nlm\.nih\.gov/articles/(PMC\d+)", identifier, re.IGNORECASE)
    if pmc_url:
        return {"type": "pmcid", "value": pmc_url.group(1).upper()}

    # Extract DOI from nature/other URLs
    doi_url = re.search(r"doi\.org/(10\.\S+)", identifier)
    if doi_url:
        return {"type": "doi", "value": doi_url.group(1).rstrip("/")}

    # Extract DOI from nature URLs like nature.com/articles/s41586-024-07259-6
    nature_url = re.search(r"nature\.com/articles/(s\d+-\d+-\d+-\w+)", identifier)
    if nature_url:
        return {"type": "doi", "value": f"10.1038/{nature_url.group(1)}"}

    # Direct PMID
    if PMID_PATTERN.match(identifier):
        return {"type": "pmid", "value": identifier}

    # Direct PMCID
    if PMCID_PATTERN.match(identifier):
        return {"type": "pmcid", "value": identifier.upper()}

    # Direct DOI
    if DOI_PATTERN.match(identifier):
        return {"type": "doi", "value": identifier}

    return {"type": "unknown", "value": identifier}


def load_identifiers_from_file(filepath: str) -> list[dict]:
    """Load paper identifiers from a text file (one per line).

    Returns list of parsed identifier dicts.
    """
    path = Path(filepath)
    if not path.exists():
        console.print(f"[red]File not found: {filepath}[/red]")
        return []

    identifiers = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            parsed = parse_identifier(line)
            if parsed["type"] != "unknown":
                identifiers.append(parsed)
            else:
                console.print(f"[yellow]Skipping unrecognized identifier: {line}[/yellow]")

    console.print(f"Loaded {len(identifiers)} identifiers from {filepath}")
    return identifiers


def resolve_to_pmid(identifier: dict) -> str | None:
    """Resolve any identifier type to a PMID using NCBI E-utilities.

    For DOIs and PMCIDs, queries NCBI to find the corresponding PMID.
    """
    from Bio import Entrez
    import config

    Entrez.email = config.NCBI_EMAIL

    id_type = identifier["type"]
    value = identifier["value"]

    if id_type == "pmid":
        return value

    if id_type == "pmcid":
        try:
            handle = Entrez.esearch(db="pubmed", term=f"{value}[pmcid]")
            result = Entrez.read(handle)
            handle.close()
            ids = result.get("IdList", [])
            return ids[0] if ids else None
        except Exception as e:
            console.print(f"[red]Error resolving PMCID {value}: {e}[/red]")
            return None

    if id_type == "doi":
        try:
            handle = Entrez.esearch(db="pubmed", term=f"{value}[doi]")
            result = Entrez.read(handle)
            handle.close()
            ids = result.get("IdList", [])
            return ids[0] if ids else None
        except Exception as e:
            console.print(f"[red]Error resolving DOI {value}: {e}[/red]")
            return None

    return None
