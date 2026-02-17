"""Manual PMID/DOI input handling."""
import re
import logging
from pathlib import Path
from urllib.parse import urlparse, urljoin
from rich.console import Console

console = Console()
logger = logging.getLogger(__name__)

# Patterns
PMID_PATTERN = re.compile(r"^\d{6,9}$")
PMCID_PATTERN = re.compile(r"^PMC\d+$", re.IGNORECASE)
DOI_PATTERN = re.compile(r"^10\.\d{4,}/\S+$")

# Publisher URL patterns that embed a DOI in the path
# Each tuple: (domain_regex, path_regex_with_doi_capture)
DOI_URL_PATTERNS = [
    # doi.org redirect
    (r"doi\.org", r"/(10\.\d{4,}/\S+)"),
    # Nature: nature.com/articles/s41586-024-07259-6
    (r"nature\.com", r"/articles/(s\d+-\d+-\d+-\w+)"),
    # Wiley: onlinelibrary.wiley.com/doi/10.xxxx/...
    (r"wiley\.com", r"/doi/(10\.\d{4,}/[^\s?#]+)"),
    # Science/AAAS: science.org/doi/10.xxxx/...
    (r"science\.org", r"/doi/(10\.\d{4,}/[^\s?#]+)"),
    # Cell Press: cell.com/*/fulltext/S0092-8674(22)01234-5 → use page fetch
    # ACS: pubs.acs.org/doi/10.xxxx/...
    (r"acs\.org", r"/doi/(10\.\d{4,}/[^\s?#]+)"),
    # bioRxiv/medRxiv: biorxiv.org/content/10.xxxx/...
    (r"biorxiv\.org|medrxiv\.org", r"/content/(10\.\d{4,}/[^\s?#v]+)"),
    # Oxford Academic: academic.oup.com/*/article-doi/10.xxxx/...
    (r"oup\.com", r"/(10\.\d{4,}/[^\s?#]+)"),
    # Springer/BMC: link.springer.com/article/10.xxxx/...
    (r"springer\.com|biomedcentral\.com", r"/article/(10\.\d{4,}/[^\s?#]+)"),
    # Frontiers: frontiersin.org/articles/10.xxxx/...
    (r"frontiersin\.org", r"/articles/(10\.\d{4,}/[^\s?#]+)"),
    # Taylor & Francis: tandfonline.com/doi/10.xxxx/...
    (r"tandfonline\.com", r"/doi/(?:abs|full)?/?(10\.\d{4,}/[^\s?#]+)"),
    # PLOS: journals.plos.org/plosone/article?id=10.xxxx/...
    (r"plos\.org", r"[?&]id=(10\.\d{4,}/[^\s&#]+)"),
    # PNAS: pnas.org/doi/10.xxxx/...
    (r"pnas\.org", r"/doi/(10\.\d{4,}/[^\s?#]+)"),
    # Elsevier/ScienceDirect: sciencedirect.com/science/article/pii/... → use page fetch
    # Generic catch-all: any URL with /10.xxxx/ in the path
    (r".", r"/(10\.\d{4,}/[^\s?#]+)"),
]


def _extract_doi_from_url(url: str) -> str | None:
    """Try to extract a DOI from a publisher URL using known patterns."""
    for domain_re, path_re in DOI_URL_PATTERNS:
        if re.search(domain_re, url, re.IGNORECASE):
            m = re.search(path_re, url)
            if m:
                doi = m.group(1).rstrip("/.,;")
                # Validate it looks like a real DOI
                if re.match(r"^10\.\d{4,}/", doi):
                    return doi
    return None


def _fetch_doi_from_page(url: str) -> str | None:
    """Fetch a URL and extract DOI from HTML meta tags.

    Looks for common citation meta tags:
      <meta name="citation_doi" content="10.xxxx/...">
      <meta name="DC.identifier" content="doi:10.xxxx/...">
      <meta name="doi" content="10.xxxx/...">
    """
    try:
        import requests
        resp = requests.get(
            url,
            headers={"User-Agent": "pegrna-db/1.0 (paper resolver)"},
            timeout=15,
            allow_redirects=True,
        )
        resp.raise_for_status()
        html = resp.text[:50_000]  # Only scan first 50KB

        # Try meta tags
        for pattern in [
            r'<meta\s+name=["\'](?:citation_doi|doi|DC\.identifier)["\']\s+content=["\'](?:doi:)?(10\.\d{4,}/[^"\']+)["\']',
            r'<meta\s+content=["\'](?:doi:)?(10\.\d{4,}/[^"\']+)["\']\s+name=["\'](?:citation_doi|doi|DC\.identifier)["\']',
            r'doi\.org/(10\.\d{4,}/[^"\'<\s]+)',
        ]:
            m = re.search(pattern, html, re.IGNORECASE)
            if m:
                return m.group(1).rstrip("/.,;")
    except Exception as e:
        logger.debug(f"Failed to fetch DOI from {url}: {e}")
    return None


def _resolve_doi_via_crossref(url: str) -> str | None:
    """Try to resolve a publisher URL to a DOI via the CrossRef API.

    Extracts ISSN / article identifiers from the URL path and queries
    CrossRef to find the corresponding DOI.
    """
    try:
        import requests

        clean = re.sub(r"[?#].*$", "", url)
        parsed = urlparse(clean)
        path_parts = [p for p in parsed.path.strip("/").split("/") if p]

        # MDPI pattern: mdpi.com/{ISSN}/{volume}/{issue}/{article}
        # e.g. mdpi.com/2073-4425/13/12/2348
        if "mdpi.com" in parsed.netloc and len(path_parts) >= 4:
            issn = path_parts[0]
            article_num = path_parts[-1]
            if re.match(r"\d{4}-\d{3}[\dX]", issn):
                resp = requests.get(
                    f"https://api.crossref.org/journals/{issn}/works",
                    params={"query": article_num, "rows": 3, "select": "DOI,title"},
                    headers={"User-Agent": "pegrna-db/1.0 (mailto:joey.riepsaame@path.ox.ac.uk)"},
                    timeout=15,
                )
                if resp.status_code == 200:
                    items = resp.json().get("message", {}).get("items", [])
                    # Find the item whose DOI ends with the article number
                    for item in items:
                        doi = item.get("DOI", "")
                        if doi.endswith(article_num):
                            return doi
                    # Fallback: return first result if only one
                    if len(items) == 1:
                        return items[0].get("DOI")

        # Generic: extract any path components that look like identifiers
        # and search CrossRef
        query_parts = [p for p in path_parts if len(p) > 3 and not p.startswith("article")]
        if query_parts:
            query = " ".join(query_parts[-3:])  # Use last 3 path components
            resp = requests.get(
                "https://api.crossref.org/works",
                params={"query": query, "rows": 1, "select": "DOI,title"},
                headers={"User-Agent": "pegrna-db/1.0 (mailto:joey.riepsaame@path.ox.ac.uk)"},
                timeout=15,
            )
            if resp.status_code == 200:
                items = resp.json().get("message", {}).get("items", [])
                if items:
                    return items[0].get("DOI")

    except Exception as e:
        logger.debug(f"CrossRef resolution failed for {url}: {e}")
    return None


def parse_identifier(identifier: str) -> dict:
    """Parse a paper identifier (PMID, PMCID, DOI, or URL).

    Supports:
      - Direct PMIDs, PMCIDs, DOIs
      - PubMed and PMC URLs
      - Publisher URLs (Nature, MDPI, Cell, Wiley, Science, bioRxiv, etc.)
      - Any URL with an embedded DOI or citation_doi meta tag

    Returns dict with keys: type, value (normalized).
    """
    identifier = identifier.strip()

    # Strip tracking parameters for cleaner matching
    clean = re.sub(r'[?#].*$', '', identifier)

    # Extract PMID from PubMed URL
    pmid_url = re.search(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)", identifier)
    if pmid_url:
        return {"type": "pmid", "value": pmid_url.group(1)}

    # Extract PMCID from PMC URL
    pmc_url = re.search(r"pmc\.ncbi\.nlm\.nih\.gov/articles/(PMC\d+)", identifier, re.IGNORECASE)
    if pmc_url:
        return {"type": "pmcid", "value": pmc_url.group(1).upper()}

    # Direct PMID
    if PMID_PATTERN.match(identifier):
        return {"type": "pmid", "value": identifier}

    # Direct PMCID
    if PMCID_PATTERN.match(identifier):
        return {"type": "pmcid", "value": identifier.upper()}

    # Direct DOI
    if DOI_PATTERN.match(identifier):
        return {"type": "doi", "value": identifier}

    # Try to extract DOI from publisher URLs
    if identifier.startswith(("http://", "https://")):
        # Nature special case (DOI prefix mapping)
        nature_url = re.search(r"nature\.com/articles/(s\d+-\d+-\d+-\w+)", clean)
        if nature_url:
            return {"type": "doi", "value": f"10.1038/{nature_url.group(1)}"}

        # Try extracting DOI from URL path (cleaned URL without query params)
        doi = _extract_doi_from_url(clean)
        if doi:
            return {"type": "doi", "value": doi}

        # Also try the full URL (for PLOS ?id=10.xxxx/... style)
        if clean != identifier:
            doi = _extract_doi_from_url(identifier)
            if doi:
                return {"type": "doi", "value": doi}

        # Fallback: fetch page / CrossRef API
        return {"type": "url", "value": identifier}

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
    For URLs, first extracts DOI from the page, then resolves to PMID.
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

    if id_type == "url":
        # Strategy 1: Fetch page and extract DOI from meta tags
        doi = _fetch_doi_from_page(value)
        if doi:
            console.print(f"[green]Extracted DOI from page: {doi}[/green]")
            return resolve_to_pmid({"type": "doi", "value": doi})

        # Strategy 2: Resolve via CrossRef API
        doi = _resolve_doi_via_crossref(value)
        if doi:
            console.print(f"[green]Resolved DOI via CrossRef: {doi}[/green]")
            return resolve_to_pmid({"type": "doi", "value": doi})

        console.print(f"[red]Could not extract DOI from URL: {value}[/red]")
        return None

    return None
