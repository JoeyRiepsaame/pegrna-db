"""PubMed search via Biopython Entrez."""
import time
from typing import Optional
from Bio import Entrez
from rich.console import Console

import config

console = Console()

Entrez.email = config.NCBI_EMAIL
Entrez.tool = config.NCBI_TOOL


def search_pubmed(
    query: Optional[str] = None,
    max_results: int = config.PUBMED_MAX_RESULTS,
) -> list[str]:
    """Search PubMed and return list of PMIDs.

    Args:
        query: Custom search query, or None to use default prime editing terms.
        max_results: Maximum number of results to return.

    Returns:
        List of PMID strings.
    """
    if query is None:
        query = config.PUBMED_SEARCH_TERMS[0]

    console.print(f"[bold]Searching PubMed:[/bold] {query}")

    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=max_results,
        sort="relevance",
    )
    results = Entrez.read(handle)
    handle.close()

    pmids = results.get("IdList", [])
    total = results.get("Count", "0")
    console.print(f"Found [green]{total}[/green] results, returning {len(pmids)}")

    return pmids


def fetch_paper_metadata(pmids: list[str]) -> list[dict]:
    """Fetch metadata for a list of PMIDs.

    Returns list of dicts with keys: pmid, title, authors, journal, year, abstract, doi, pmcid.
    """
    if not pmids:
        return []

    papers = []
    # Process in batches of 50 to respect rate limits
    batch_size = 50
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i : i + batch_size]
        console.print(f"Fetching metadata for PMIDs {i+1}-{i+len(batch)} of {len(pmids)}...")

        handle = Entrez.efetch(
            db="pubmed",
            id=",".join(batch),
            rettype="xml",
            retmode="xml",
        )
        records = Entrez.read(handle)
        handle.close()

        for article in records.get("PubmedArticle", []):
            paper = _parse_pubmed_article(article)
            if paper:
                papers.append(paper)

        if i + batch_size < len(pmids):
            time.sleep(0.4)  # Rate limit: ~3 req/sec without API key

    return papers


def _parse_pubmed_article(article: dict) -> Optional[dict]:
    """Parse a PubMed article record into a flat dict."""
    try:
        medline = article.get("MedlineCitation", {})
        article_data = medline.get("Article", {})
        pmid = str(medline.get("PMID", ""))

        # Title
        title = str(article_data.get("ArticleTitle", ""))

        # Authors
        author_list = article_data.get("AuthorList", [])
        authors = []
        for a in author_list:
            last = a.get("LastName", "")
            first = a.get("ForeName", "")
            if last:
                authors.append(f"{last} {first}".strip())
        authors_str = ", ".join(authors)

        # Journal
        journal_info = article_data.get("Journal", {})
        journal = str(journal_info.get("Title", ""))

        # Year
        pub_date = journal_info.get("JournalIssue", {}).get("PubDate", {})
        year = pub_date.get("Year")
        if year:
            year = int(year)
        else:
            # Try MedlineDate
            medline_date = pub_date.get("MedlineDate", "")
            if medline_date and medline_date[:4].isdigit():
                year = int(medline_date[:4])

        # Abstract
        abstract_parts = article_data.get("Abstract", {}).get("AbstractText", [])
        abstract = " ".join(str(part) for part in abstract_parts)

        # DOI
        doi = None
        id_list = article_data.get("ELocationID", [])
        for eid in id_list:
            if eid.attributes.get("EIdType") == "doi":
                doi = str(eid)
                break

        # PMC ID
        pmcid = None
        pmc_data = article.get("PubmedData", {})
        for aid in pmc_data.get("ArticleIdList", []):
            if aid.attributes.get("IdType") == "pmc":
                pmcid = str(aid)
                break

        return {
            "pmid": pmid,
            "title": title,
            "authors": authors_str,
            "journal": journal,
            "year": year,
            "abstract": abstract,
            "doi": doi,
            "pmcid": pmcid,
        }
    except Exception as e:
        console.print(f"[red]Error parsing article: {e}[/red]")
        return None


def pmid_to_pmcid(pmid: str) -> Optional[str]:
    """Convert PMID to PMCID using ID converter."""
    try:
        handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
        result = Entrez.read(handle)
        handle.close()

        for linkset in result:
            for link_db in linkset.get("LinkSetDb", []):
                if link_db.get("DbTo") == "pmc":
                    links = link_db.get("Link", [])
                    if links:
                        return f"PMC{links[0]['Id']}"
    except Exception:
        pass
    return None
