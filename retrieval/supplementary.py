"""Download and parse supplementary materials from PMC."""
import io
import tarfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional

import pandas as pd
import requests
from rich.console import Console

import config

console = Console()


def list_supplementary_files(pmcid: str) -> list[dict]:
    """List available supplementary files for a PMC article.

    Uses PMC OA API to find the tar.gz package, then lists its contents.
    Returns list of dicts with keys: filename, type.
    """
    # Use OA API to get the FTP link for the package
    oa_url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmcid}"
    try:
        resp = requests.get(oa_url, timeout=15)
        if resp.status_code != 200:
            console.print(f"[yellow]OA API returned {resp.status_code} for {pmcid}[/yellow]")
            return []

        root = ET.fromstring(resp.text)
        records = root.findall(".//record")
        if not records:
            console.print(f"[yellow]No OA package for {pmcid}[/yellow]")
            return []

        # Get the tgz link
        tgz_link = None
        for record in records:
            for link in record.findall("link"):
                if link.get("format") == "tgz":
                    tgz_link = link.get("href")
                    break

        if not tgz_link:
            console.print(f"[yellow]No tgz package link for {pmcid}[/yellow]")
            return []

        # Convert ftp to https
        tgz_link = tgz_link.replace("ftp://ftp.ncbi.nlm.nih.gov/", "https://ftp.ncbi.nlm.nih.gov/")

        console.print(f"Downloading package for {pmcid}...")
        pkg_resp = requests.get(tgz_link, timeout=120)
        if pkg_resp.status_code != 200:
            console.print(f"[yellow]Failed to download package (HTTP {pkg_resp.status_code})[/yellow]")
            return []

        # Extract and list files
        files = []
        save_dir = config.RAW_PAPERS_DIR / pmcid
        save_dir.mkdir(parents=True, exist_ok=True)

        with tarfile.open(fileobj=io.BytesIO(pkg_resp.content), mode="r:gz") as tar:
            for member in tar.getmembers():
                if member.isfile():
                    fname = Path(member.name).name
                    ftype = _classify_file_type(fname)
                    files.append({"filename": fname, "type": ftype})

                    # Extract tabular files
                    if ftype in ("excel", "csv", "tsv"):
                        tar.extract(member, path=str(save_dir))
                        # Move from nested dir to save_dir
                        extracted = save_dir / member.name
                        target = save_dir / fname
                        if extracted != target and extracted.exists():
                            extracted.rename(target)

        console.print(f"Found {len(files)} files in {pmcid} package")
        return files

    except Exception as e:
        console.print(f"[red]Error listing supplementary files: {e}[/red]")
        return []


def download_supplementary_file(
    pmcid: str, filename: str, save_dir: Optional[Path] = None,
) -> Optional[Path]:
    """Get path to an already-extracted supplementary file.

    Files are extracted during list_supplementary_files().
    Falls back to direct download from PMC bin directory.
    """
    if save_dir is None:
        save_dir = config.RAW_PAPERS_DIR / pmcid

    filepath = save_dir / filename
    if filepath.exists():
        console.print(f"[green]Found: {filename}[/green]")
        return filepath

    # Fallback: try direct download from PMC bin
    base_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/bin/"
    file_url = f"{base_url}{filename}"

    try:
        save_dir.mkdir(parents=True, exist_ok=True)
        resp = requests.get(file_url, timeout=60)
        if resp.status_code == 200:
            filepath.write_bytes(resp.content)
            console.print(f"[green]Downloaded: {filename}[/green]")
            return filepath
        else:
            console.print(f"[yellow]Failed to download {filename} (HTTP {resp.status_code})[/yellow]")
    except Exception as e:
        console.print(f"[red]Error downloading {filename}: {e}[/red]")

    return None


def parse_supplementary_tables(filepath: Path) -> list[pd.DataFrame]:
    """Parse a supplementary file into pandas DataFrames.

    Supports: .xlsx, .xls, .csv, .tsv
    Auto-detects header rows for Excel files where headers aren't on row 0.
    Returns list of DataFrames (one per sheet for Excel files).
    """
    suffix = filepath.suffix.lower()
    dfs = []

    try:
        if suffix in (".xlsx", ".xls"):
            xl = pd.ExcelFile(filepath)
            for sheet_name in xl.sheet_names:
                try:
                    # First read without header to detect the real header row
                    raw = pd.read_excel(xl, sheet_name=sheet_name, header=None)
                    if raw.empty:
                        continue

                    header_row = _find_header_row(raw)
                    df = pd.read_excel(xl, sheet_name=sheet_name, header=header_row)

                    # Drop rows that are all NaN
                    df = df.dropna(how="all")
                    if not df.empty:
                        df.attrs["source_sheet"] = sheet_name
                        df.attrs["source_file"] = filepath.name
                        dfs.append(df)
                except Exception:
                    continue

        elif suffix == ".csv":
            df = pd.read_csv(filepath)
            if not df.empty:
                df.attrs["source_file"] = filepath.name
                dfs.append(df)

        elif suffix in (".tsv", ".txt"):
            df = pd.read_csv(filepath, sep="\t")
            if not df.empty:
                df.attrs["source_file"] = filepath.name
                dfs.append(df)

    except Exception as e:
        console.print(f"[red]Error parsing {filepath.name}: {e}[/red]")

    console.print(
        f"Parsed {len(dfs)} table(s) from {filepath.name} "
        f"({sum(len(df) for df in dfs)} total rows)"
    )
    return dfs


def _find_header_row(df: pd.DataFrame, max_scan: int = 20) -> int:
    """Find the likely header row in a DataFrame.

    Looks for rows with multiple short string values that look like column names.
    Checks for known pegRNA-related terms.
    """
    known_headers = {
        "spacer", "pbs", "rtt", "template", "scaffold", "guide", "sequence",
        "efficiency", "gene", "target", "locus", "motif", "linker", "name",
        "edit", "cell", "type", "primer", "amplicon", "full length",
    }

    best_row = 0
    best_score = 0

    for i in range(min(max_scan, len(df))):
        row_values = [str(v).lower().strip() for v in df.iloc[i] if pd.notna(v)]
        if not row_values:
            continue

        # Score: how many values match known headers
        score = sum(
            1 for v in row_values
            if any(h in v for h in known_headers)
        )
        # Bonus for having multiple non-empty short values (header-like)
        short_values = [v for v in row_values if 0 < len(v) < 50]
        score += len(short_values) * 0.1

        if score > best_score:
            best_score = score
            best_row = i

    return best_row


def find_pegrna_tables(dfs: list[pd.DataFrame]) -> list[pd.DataFrame]:
    """Filter DataFrames to those likely containing pegRNA data.

    Checks column names against known pegRNA-related patterns.
    """
    from config import (
        SPACER_PATTERNS, PBS_PATTERNS, RTT_PATTERNS,
        EFFICIENCY_PATTERNS, GENE_PATTERNS,
    )

    all_patterns = (
        SPACER_PATTERNS + PBS_PATTERNS + RTT_PATTERNS
        + EFFICIENCY_PATTERNS + GENE_PATTERNS
    )

    relevant = []
    for df in dfs:
        cols_lower = [str(c).lower().strip() for c in df.columns]
        matches = sum(
            1 for col in cols_lower
            if any(pat in col for pat in all_patterns)
        )
        if matches >= 2:
            relevant.append(df)
            sheet = df.attrs.get("source_sheet", "")
            file = df.attrs.get("source_file", "")
            console.print(
                f"[green]Found pegRNA table: {file}"
                f"{f' / {sheet}' if sheet else ''} "
                f"({len(df)} rows, {matches} matching columns)[/green]"
            )

    return relevant


def _classify_file_type(filename: str) -> str:
    """Classify supplementary file type."""
    suffix = Path(filename).suffix.lower()
    if suffix in (".xlsx", ".xls"):
        return "excel"
    elif suffix == ".csv":
        return "csv"
    elif suffix in (".tsv", ".txt"):
        return "tsv"
    elif suffix == ".pdf":
        return "pdf"
    elif suffix in (".zip", ".gz", ".tar"):
        return "archive"
    else:
        return "other"
