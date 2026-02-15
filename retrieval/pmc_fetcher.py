"""Fetch full text from PMC via BioC API."""
import json
import re
import time
from html.parser import HTMLParser
from typing import Optional
import requests
from rich.console import Console

import config

console = Console()


class _HTMLTextExtractor(HTMLParser):
    """Simple HTML tag stripper."""
    def __init__(self):
        super().__init__()
        self._texts = []
        self._skip = False

    def handle_starttag(self, tag, attrs):
        if tag in ("script", "style", "nav", "header", "footer"):
            self._skip = True

    def handle_endtag(self, tag):
        if tag in ("script", "style", "nav", "header", "footer"):
            self._skip = False

    def handle_data(self, data):
        if not self._skip:
            self._texts.append(data)

    def get_text(self):
        return " ".join(self._texts)


def fetch_full_text(pmcid: str, max_retries: int = 3) -> Optional[dict]:
    """Fetch full text of a PMC article in BioC JSON format.

    Args:
        pmcid: PMC ID (e.g., "PMC8930418").
        max_retries: Number of retry attempts on transient errors.

    Returns:
        Parsed BioC JSON dict, or None if unavailable.
    """
    pmc_num = pmcid.replace("PMC", "")
    url = f"{config.PMC_BIOC_BASE}/{pmcid}/unicode"

    console.print(f"Fetching full text for {pmcid}...")
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                if isinstance(data, list) and data:
                    data = data[0]
                console.print(f"[green]Got full text for {pmcid}[/green]")
                return data
            elif resp.status_code in (429, 503):
                wait = 2 ** attempt * 5
                console.print(f"[yellow]Rate limited, waiting {wait}s...[/yellow]")
                time.sleep(wait)
                continue
            else:
                console.print(
                    f"[yellow]PMC BioC returned {resp.status_code} for {pmcid}[/yellow]"
                )
                return None
        except json.JSONDecodeError:
            console.print(f"[yellow]BioC JSON decode failed for {pmcid}[/yellow]")
            return None
        except Exception as e:
            console.print(f"[red]Error fetching {pmcid}: {e}[/red]")
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)
    return None


def fetch_full_text_html(pmcid: str) -> Optional[str]:
    """Fallback: fetch PMC article as HTML and extract plain text.

    For papers where BioC API returns empty or error.
    """
    url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
    console.print(f"Trying HTML fallback for {pmcid}...")
    try:
        resp = requests.get(
            url, timeout=30,
            headers={"User-Agent": "pegrna-db/1.0 (research)"},
        )
        if resp.status_code == 200 and len(resp.content) > 5000:
            parser = _HTMLTextExtractor()
            parser.feed(resp.text)
            text = parser.get_text()
            if len(text) > 500:
                console.print(f"[green]Got HTML text for {pmcid} ({len(text)} chars)[/green]")
                return text
    except Exception as e:
        console.print(f"[red]HTML fallback failed for {pmcid}: {e}[/red]")
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
