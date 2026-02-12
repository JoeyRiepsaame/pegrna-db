"""Unpaywall API for open access PDF fallback."""
from typing import Optional
import requests
from rich.console import Console

import config

console = Console()


def get_open_access_url(doi: str) -> Optional[str]:
    """Look up OA PDF URL via Unpaywall API.

    Args:
        doi: DOI of the paper.

    Returns:
        URL string to the best OA version, or None.
    """
    if not doi or not config.UNPAYWALL_EMAIL:
        return None

    url = f"{config.UNPAYWALL_BASE}/{doi}"
    params = {"email": config.UNPAYWALL_EMAIL}

    try:
        resp = requests.get(url, params=params, timeout=15)
        if resp.status_code == 200:
            data = resp.json()
            best = data.get("best_oa_location")
            if best:
                pdf_url = best.get("url_for_pdf") or best.get("url")
                if pdf_url:
                    console.print(f"[green]Found OA link via Unpaywall: {pdf_url}[/green]")
                    return pdf_url

        console.print(f"[yellow]No OA version found for DOI {doi}[/yellow]")
    except Exception as e:
        console.print(f"[red]Unpaywall error: {e}[/red]")

    return None
