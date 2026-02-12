"""Fetch full text from PMC via BioC API."""
import json
import time
from typing import Optional
import requests
from rich.console import Console

import config

console = Console()


def fetch_full_text(pmcid: str) -> Optional[dict]:
    """Fetch full text of a PMC article in BioC JSON format.

    Args:
        pmcid: PMC ID (e.g., "PMC8930418").

    Returns:
        Parsed BioC JSON dict, or None if unavailable.
    """
    # Strip 'PMC' prefix if present for the API
    pmc_num = pmcid.replace("PMC", "")
    url = f"{config.PMC_BIOC_BASE}/{pmcid}/unicode"

    console.print(f"Fetching full text for {pmcid}...")
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            # BioC API returns a list containing one dict
            if isinstance(data, list) and data:
                data = data[0]
            console.print(f"[green]Got full text for {pmcid}[/green]")
            return data
        else:
            console.print(
                f"[yellow]PMC BioC returned {resp.status_code} for {pmcid}[/yellow]"
            )
            return None
    except Exception as e:
        console.print(f"[red]Error fetching {pmcid}: {e}[/red]")
        return None


def extract_text_from_bioc(bioc_data: dict) -> str:
    """Extract plain text from BioC JSON structure.

    Returns concatenated text from all passages.
    """
    texts = []
    for doc in bioc_data.get("documents", []):
        for passage in doc.get("passages", []):
            text = passage.get("text", "")
            section = passage.get("infons", {}).get("section_type", "")
            if text:
                if section:
                    texts.append(f"\n[{section}]\n{text}")
                else:
                    texts.append(text)
    return "\n".join(texts)


def extract_tables_from_bioc(bioc_data: dict) -> list[dict]:
    """Extract tables from BioC JSON data.

    Returns list of dicts with keys: table_id, caption, content.
    """
    tables = []
    for doc in bioc_data.get("documents", []):
        for passage in doc.get("passages", []):
            infons = passage.get("infons", {})
            section_type = infons.get("section_type", "").lower()
            passage_type = infons.get("type", "").lower()

            if "table" in section_type or "table" in passage_type:
                table = {
                    "table_id": infons.get("id", ""),
                    "caption": "",
                    "content": passage.get("text", ""),
                }
                if "caption" in passage_type or "title" in passage_type:
                    table["caption"] = passage.get("text", "")
                tables.append(table)

    return tables


def fetch_by_pmid(pmid: str) -> Optional[dict]:
    """Try to fetch full text using PMID (converts to PMCID first)."""
    from discovery.pubmed_search import pmid_to_pmcid

    pmcid = pmid_to_pmcid(pmid)
    if pmcid:
        return fetch_full_text(pmcid)

    console.print(f"[yellow]No PMC version found for PMID {pmid}[/yellow]")
    return None
