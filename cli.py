"""CLI entry point for pegrna-db."""
import sys
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

# Ensure project root is on path
sys.path.insert(0, str(Path(__file__).parent))

import config
from database.models import init_db
from database.operations import (
    get_or_create_paper, bulk_add_entries, update_paper_status,
    search_entries, get_stats, bulk_update_sequences_for_paper,
)

console = Console()
app = typer.Typer(
    name="pegrna-db",
    help="Automated (e)pegRNA database from scientific literature.",
)


def get_session():
    """Get a database session."""
    Session = init_db(str(config.DATABASE_PATH))
    return Session()


@app.command()
def search(
    query: Optional[str] = typer.Argument(
        None, help="Custom PubMed query. Uses default prime editing query if not provided."
    ),
    max_results: int = typer.Option(50, "--max", "-m", help="Maximum papers to fetch"),
    process: bool = typer.Option(True, "--process/--no-process", help="Process papers after discovery"),
):
    """Search PubMed for prime editing papers and optionally process them."""
    from discovery.pubmed_search import search_pubmed, fetch_paper_metadata

    pmids = search_pubmed(query=query, max_results=max_results)
    if not pmids:
        console.print("[yellow]No papers found[/yellow]")
        raise typer.Exit()

    papers = fetch_paper_metadata(pmids)
    session = get_session()

    added = 0
    for paper_data in papers:
        paper, created = get_or_create_paper(session, **paper_data)
        if created:
            added += 1
    session.commit()

    console.print(f"[green]Added {added} new papers[/green] ({len(papers) - added} already in database)")

    if process and added > 0:
        console.print("\n[bold]Processing new papers...[/bold]")
        _process_pending_papers(session)

    session.close()


@app.command()
def add(
    identifiers: list[str] = typer.Argument(help="PMIDs, DOIs, or URLs to add"),
    process: bool = typer.Option(True, "--process/--no-process", help="Process after adding"),
):
    """Add specific papers by PMID, DOI, or URL."""
    from discovery.manual_input import parse_identifier, resolve_to_pmid
    from discovery.pubmed_search import fetch_paper_metadata

    pmids = []
    for ident in identifiers:
        parsed = parse_identifier(ident)
        console.print(f"Parsed: {parsed['type']} = {parsed['value']}")

        if parsed["type"] == "pmid":
            pmids.append(parsed["value"])
        elif parsed["type"] in ("doi", "pmcid"):
            pmid = resolve_to_pmid(parsed)
            if pmid:
                pmids.append(pmid)
            else:
                console.print(f"[red]Could not resolve {ident} to PMID[/red]")
        else:
            console.print(f"[red]Unrecognized identifier: {ident}[/red]")

    if not pmids:
        console.print("[red]No valid identifiers found[/red]")
        raise typer.Exit(1)

    papers = fetch_paper_metadata(pmids)
    session = get_session()

    added = 0
    for paper_data in papers:
        paper, created = get_or_create_paper(session, **paper_data)
        if created:
            added += 1
            console.print(f"  [green]+[/green] {paper.title}")
        else:
            console.print(f"  [yellow]=[/yellow] Already exists: {paper.title}")
    session.commit()

    if process:
        _process_pending_papers(session)

    session.close()


@app.command()
def extract(
    paper_id: Optional[int] = typer.Argument(None, help="Specific paper ID to re-extract"),
    method: str = typer.Option("auto", help="Extraction method: auto, rule, llm"),
):
    """Run extraction on pending papers or a specific paper."""
    session = get_session()

    if paper_id:
        from database.models import Paper
        paper = session.query(Paper).get(paper_id)
        if not paper:
            console.print(f"[red]Paper ID {paper_id} not found[/red]")
            raise typer.Exit(1)
        _process_single_paper(session, paper, force_method=method)
    else:
        _process_pending_papers(session, force_method=method)

    session.close()


@app.command()
def stats():
    """Show database statistics."""
    session = get_session()
    s = get_stats(session)

    table = Table(title="pegRNA Database Statistics")
    table.add_column("Metric", style="bold")
    table.add_column("Value", justify="right")

    table.add_row("Total papers", str(s["total_papers"]))
    table.add_row("Processed papers", str(s["completed_papers"]))
    table.add_row("Total pegRNA entries", str(s["total_entries"]))
    table.add_row("Validated entries", str(s["validated_entries"]))
    table.add_row("Unique genes", str(s["unique_genes"]))
    table.add_row("Unique cell types", str(s["unique_cell_types"]))
    table.add_row("Unique organisms", str(s["unique_organisms"]))

    console.print(table)
    session.close()


@app.command()
def export(
    format: str = typer.Option("csv", "--format", "-f", help="Output format: csv, json"),
    output: str = typer.Option("pegrna_export", "--output", "-o", help="Output filename (without extension)"),
    validated_only: bool = typer.Option(False, "--validated", help="Export only validated entries"),
):
    """Export database to CSV or JSON."""
    import pandas as pd
    from database.models import PegRNAEntry, Paper

    session = get_session()
    query = session.query(PegRNAEntry).join(Paper)

    if validated_only:
        query = query.filter(PegRNAEntry.validated.is_(True))

    entries = query.all()
    if not entries:
        console.print("[yellow]No entries to export[/yellow]")
        raise typer.Exit()

    rows = []
    for e in entries:
        row = {
            "entry_name": e.entry_name,
            "pegrna_type": e.pegrna_type,
            "spacer_sequence": e.spacer_sequence,
            "pbs_sequence": e.pbs_sequence,
            "pbs_length": e.pbs_length,
            "rtt_sequence": e.rtt_sequence,
            "rtt_length": e.rtt_length,
            "three_prime_extension": e.three_prime_extension,
            "full_sequence": e.full_sequence,
            "nicking_sgrna_seq": e.nicking_sgrna_seq,
            "target_gene": e.target_gene,
            "target_locus": e.target_locus,
            "target_organism": e.target_organism,
            "edit_type": e.edit_type,
            "edit_description": e.edit_description,
            "intended_mutation": e.intended_mutation,
            "prime_editor": e.prime_editor,
            "cell_type": e.cell_type,
            "delivery_method": e.delivery_method,
            "editing_efficiency": e.editing_efficiency,
            "product_purity": e.product_purity,
            "indel_frequency": e.indel_frequency,
            "confidence_score": e.confidence_score,
            "validated": e.validated,
            "paper_pmid": e.paper.pmid,
            "paper_title": e.paper.title,
            "paper_doi": e.paper.doi,
            "paper_year": e.paper.year,
        }
        rows.append(row)

    df = pd.DataFrame(rows)

    if format == "csv":
        filepath = f"{output}.csv"
        df.to_csv(filepath, index=False)
    elif format == "json":
        filepath = f"{output}.json"
        df.to_json(filepath, orient="records", indent=2)
    else:
        console.print(f"[red]Unknown format: {format}[/red]")
        raise typer.Exit(1)

    console.print(f"[green]Exported {len(df)} entries to {filepath}[/green]")
    session.close()


@app.command()
def validate():
    """Run validation on all entries and report issues."""
    from extraction.validator import validate_entry
    from database.models import PegRNAEntry

    session = get_session()
    entries = session.query(PegRNAEntry).all()

    if not entries:
        console.print("[yellow]No entries in database[/yellow]")
        raise typer.Exit()

    valid_count = 0
    warning_count = 0
    invalid_count = 0

    for entry in entries:
        from extraction.schemas import PegRNAExtracted
        # Convert DB entry to schema for validation
        schema_entry = PegRNAExtracted(
            spacer_sequence=entry.spacer_sequence,
            pbs_sequence=entry.pbs_sequence,
            pbs_length=entry.pbs_length,
            rtt_sequence=entry.rtt_sequence,
            rtt_length=entry.rtt_length,
            full_sequence=entry.full_sequence,
            target_gene=entry.target_gene,
            editing_efficiency=entry.editing_efficiency,
            pegrna_type=entry.pegrna_type,
            three_prime_extension=entry.three_prime_extension,
            entry_name=entry.entry_name,
        )
        is_valid, warnings = validate_entry(schema_entry)
        if is_valid and not warnings:
            valid_count += 1
        elif is_valid:
            warning_count += 1
        else:
            invalid_count += 1
            console.print(
                f"[red]Invalid entry #{entry.id} ({entry.entry_name}): "
                f"{'; '.join(warnings)}[/red]"
            )

    console.print(
        f"\nValidation: [green]{valid_count} valid[/green], "
        f"[yellow]{warning_count} warnings[/yellow], "
        f"[red]{invalid_count} invalid[/red]"
    )
    session.close()


@app.command()
def serve(
    port: int = typer.Option(8501, help="Port for Streamlit app"),
):
    """Launch the Streamlit web app."""
    import subprocess
    webapp_path = Path(__file__).parent / "app.py"
    subprocess.run(
        ["streamlit", "run", str(webapp_path), "--server.port", str(port)],
        cwd=str(Path(__file__).parent),
    )


# --- Internal processing functions ---

def _process_pending_papers(session, force_method: str = "auto"):
    """Process all papers with pending extraction status."""
    from database.models import Paper

    pending = session.query(Paper).filter_by(extraction_status="pending").all()
    if not pending:
        console.print("[yellow]No pending papers to process[/yellow]")
        return

    console.print(f"Processing {len(pending)} pending papers...")
    for paper in pending:
        _process_single_paper(session, paper, force_method=force_method)
        session.commit()


def _process_single_paper(session, paper, force_method: str = "auto"):
    """Process a single paper through the extraction pipeline."""
    from retrieval.pmc_fetcher import fetch_full_text, extract_text_from_bioc, fetch_by_pmid, fetch_full_text_html
    from retrieval.supplementary import (
        list_supplementary_files, download_supplementary_file,
        parse_supplementary_tables, find_pegrna_tables,
    )
    from extraction.rule_based import extract_from_multiple_tables, extract_from_text
    from extraction.llm_extractor import extract_with_llm
    from extraction.validator import validate_batch

    console.print(f"\n[bold]--- Processing: {paper.title} ---[/bold]")
    console.print(f"PMID: {paper.pmid} | PMCID: {paper.pmcid} | DOI: {paper.doi}")

    all_entries = []
    paper_text = ""

    # Step 1: Try to get full text
    bioc_data = None
    if paper.pmcid:
        bioc_data = fetch_full_text(paper.pmcid)
    elif paper.pmid:
        bioc_data = fetch_by_pmid(paper.pmid)

    if bioc_data:
        paper_text = extract_text_from_bioc(bioc_data)

    # Step 1b: HTML fallback if BioC failed
    if not paper_text and paper.pmcid:
        paper_text = fetch_full_text_html(paper.pmcid) or ""

    # Step 1c: Abstract fallback when no full text available
    if not paper_text and paper.abstract:
        paper_text = f"Title: {paper.title or ''}\n\nAbstract: {paper.abstract}"
        console.print("[yellow]No full text available, using abstract for LLM extraction[/yellow]")

    # Step 2: Try supplementary materials (rule-based)
    if paper.pmcid and force_method in ("auto", "rule"):
        supp_files = list_supplementary_files(paper.pmcid)
        tabular_files = [f for f in supp_files if f["type"] in ("excel", "csv", "tsv")]

        for supp in tabular_files:
            filepath = download_supplementary_file(paper.pmcid, supp["filename"])
            if filepath:
                dfs = parse_supplementary_tables(filepath)
                pegrna_tables = find_pegrna_tables(dfs)
                if pegrna_tables:
                    entries = extract_from_multiple_tables(pegrna_tables)
                    all_entries.extend(entries)
                elif dfs and force_method in ("auto", "llm"):
                    # Tables exist but columns not recognized - try LLM
                    for df in dfs[:10]:  # Check more sheets for pegRNA data
                        table_text = df.to_string(max_rows=200)
                        from extraction.llm_extractor import extract_table_with_llm
                        entries = extract_table_with_llm(table_text, paper.title)
                        all_entries.extend(entries)

    # Step 3: Text-based extraction
    if paper_text and force_method in ("auto", "rule"):
        text_entries = extract_from_text(paper_text)
        all_entries.extend(text_entries)

    # Step 4: LLM fallback if rule-based found nothing
    if not all_entries and paper_text and force_method in ("auto", "llm"):
        console.print("[bold]Rule-based extraction found nothing, trying LLM...[/bold]")
        llm_entries = extract_with_llm(
            paper_text,
            paper_title=paper.title or "",
            paper_abstract=paper.abstract or "",
        )
        all_entries.extend(llm_entries)

    # Step 5: Validate and store
    if all_entries:
        valid_entries, invalid = validate_batch(all_entries)

        if valid_entries:
            entry_dicts = [e.to_db_dict() for e in valid_entries]
            # Remove fields that don't map to DB columns
            for d in entry_dicts:
                d.pop("raw_source_text", None)  # Keep this, it's in the model
            count = bulk_add_entries(session, paper.id, entry_dicts)
            method = "rule_based" if force_method == "rule" else "auto"
            update_paper_status(session, paper.id, "completed", method=method)
            console.print(f"[green]Stored {count} entries for {paper.pmid}[/green]")
        else:
            update_paper_status(
                session, paper.id, "review",
                notes=f"All {len(all_entries)} extracted entries failed validation",
            )
    else:
        update_paper_status(
            session, paper.id, "failed",
            notes="No pegRNA data could be extracted",
        )
        console.print(f"[red]No data extracted from {paper.pmid}[/red]")


@app.command()
def dedup():
    """Scan for duplicate papers and data quality issues."""
    from difflib import SequenceMatcher
    from database.models import Paper, PegRNAEntry
    from sqlalchemy import func, distinct

    session = get_session()
    papers = session.query(Paper).all()

    console.print("[bold]Scanning for duplicate papers...[/bold]")
    seen = []
    for paper in papers:
        if not paper.title:
            continue
        for other in seen:
            sim = SequenceMatcher(
                None, paper.title.lower(), other.title.lower()
            ).ratio()
            if sim > 0.90:
                c1 = session.query(PegRNAEntry).filter_by(paper_id=paper.id).count()
                c2 = session.query(PegRNAEntry).filter_by(paper_id=other.id).count()
                console.print(
                    f"[yellow]Potential dup (sim={sim:.2f}):[/yellow]\n"
                    f"  PMID {paper.pmid} ({c1} entries): {paper.title[:60]}\n"
                    f"  PMID {other.pmid} ({c2} entries): {other.title[:60]}"
                )
        seen.append(paper)

    console.print("\n[bold]Organism values:[/bold]")
    organisms = session.query(
        PegRNAEntry.target_organism, func.count(PegRNAEntry.id)
    ).filter(
        PegRNAEntry.target_organism.isnot(None)
    ).group_by(PegRNAEntry.target_organism).order_by(
        func.count(PegRNAEntry.id).desc()
    ).all()
    for org, count in organisms:
        console.print(f"  {org}: {count}")
    session.close()


@app.command()
def batch(
    max_papers: int = typer.Option(200, "--max", "-m"),
    skip_reviews: bool = typer.Option(True, "--skip-reviews/--include-reviews"),
    dry_run: bool = typer.Option(False, "--dry-run"),
):
    """Batch process: search PubMed, filter, and extract all new papers."""
    import time
    import json
    from discovery.pubmed_search import search_pubmed, fetch_paper_metadata
    from database.models import Paper

    REVIEW_KEYWORDS = [
        "review", "perspective", "commentary", "editorial", "meta-analysis",
        "methods only",
    ]
    EXPERIMENTAL_INDICATORS = [
        "we designed", "we tested", "we constructed", "we demonstrated",
        "editing efficiency", "we generated", "we performed",
        "we developed", "we optimized", "we applied", "we screened",
        "we achieved", "we corrected", "we established", "we validated",
        "pegRNA", "epegRNA", "editing rate", "indel", "prime edit",
        "HEK293", "U2OS", "K562", "iPSC", "fibroblast",
        "transfect", "electropor", "lentivir", "plasmid",
        "supplementary table", "supplementary data",
    ]

    session = get_session()

    # Search PubMed with multiple queries
    queries = [
        '"prime editing" AND (pegRNA OR epegRNA OR "prime editing guide RNA")',
        '"prime editing" AND (efficiency OR screen) AND (pegRNA OR guide)',
        '"prime editor" AND (PE2 OR PE3 OR PEmax OR PE5) AND pegRNA',
        '"prime editing" AND (spacer OR "guide RNA" OR "protospacer") AND efficiency',
        '"prime editing" AND (PRIDICT OR PrimeDesign OR PEGG OR pegLIT OR PE-Designer)',
        '"prime editing" AND ("high-throughput" OR "library" OR "saturation") AND screen',
    ]

    all_pmids = set()
    for q in queries:
        pmids = search_pubmed(query=q, max_results=max_papers)
        all_pmids.update(pmids)
        time.sleep(0.5)

    # Filter out existing
    existing = {p.pmid for p in session.query(Paper).all()}
    new_pmids = list(all_pmids - existing)[:max_papers]
    console.print(f"Found {len(all_pmids)} total, {len(new_pmids)} new papers to process")

    if not new_pmids:
        console.print("[green]No new papers to process[/green]")
        session.close()
        return

    # Fetch metadata
    papers_data = fetch_paper_metadata(new_pmids)
    console.print(f"Fetched metadata for {len(papers_data)} papers")

    # Filter reviews (no longer skip papers without PMCID - use abstract fallback)
    experimental = []
    skipped = 0
    no_pmcid_count = 0
    for pd in papers_data:
        if not pd.get("pmcid"):
            no_pmcid_count += 1
        if skip_reviews:
            combined = ((pd.get("title") or "") + " " + (pd.get("abstract") or "")).lower()
            is_review = any(kw in combined for kw in REVIEW_KEYWORDS)
            has_exp = any(ei in combined for ei in EXPERIMENTAL_INDICATORS)
            if is_review and not has_exp:
                skipped += 1
                continue
        experimental.append(pd)

    console.print(f"After filtering: {len(experimental)} experimental papers ({skipped} skipped, {no_pmcid_count} without PMCID)")

    if dry_run:
        for pd in experimental:
            console.print(f"  [green]+[/green] {pd.get('pmid')} - {pd.get('title', '')[:70]}")
        session.close()
        return

    # Process
    progress_path = config.DATA_DIR / "batch_progress.json"
    progress = {"processed": [], "failed": [], "entries_added": 0}

    for i, pd in enumerate(experimental):
        console.print(f"\n[bold][{i+1}/{len(experimental)}] {pd.get('title', '')[:60]}[/bold]")
        paper, created = get_or_create_paper(session, **pd)
        session.commit()

        if not created and paper.extraction_status == "completed":
            console.print(f"  [yellow]Already completed[/yellow]")
            continue

        try:
            _process_single_paper(session, paper)
            session.commit()
            entry_count = session.query(
                __import__("database.models", fromlist=["PegRNAEntry"]).PegRNAEntry
            ).filter_by(paper_id=paper.id).count()
            progress["processed"].append(pd.get("pmid"))
            progress["entries_added"] += entry_count
        except Exception as e:
            console.print(f"  [red]Error: {e}[/red]")
            progress["failed"].append(pd.get("pmid"))
            session.rollback()

        # Save progress
        with open(progress_path, "w") as f:
            json.dump(progress, f, indent=2)

    console.print(f"\n[bold green]Batch complete![/bold green]")
    console.print(f"Processed: {len(progress['processed'])}, Failed: {len(progress['failed'])}")
    console.print(f"Total new entries: {progress['entries_added']}")

    total = session.query(Paper).count()
    total_entries = session.query(
        __import__("database.models", fromlist=["PegRNAEntry"]).PegRNAEntry
    ).count()
    console.print(f"Database: {total} papers, {total_entries:,} entries")
    session.close()


@app.command(name="clinvar-update")
def clinvar_update(
    download: bool = typer.Option(True, "--download/--no-download", help="Download ClinVar data"),
    do_match: bool = typer.Option(True, "--match/--no-match", help="Run pegRNA-ClinVar matching"),
    pathogenic_only: bool = typer.Option(True, "--pathogenic-only/--all-significance",
                                         help="Only match pathogenic/likely pathogenic variants"),
):
    """Download ClinVar data, import variants, and match to pegRNA entries."""
    from database.clinvar import download_clinvar, import_clinvar_variants, match_pegrna_to_clinvar

    session = get_session()

    if download:
        filepath = download_clinvar(config.CLINVAR_DATA_DIR, config.CLINVAR_FTP_URL)
        console.print(f"\n[bold]Importing ClinVar variants...[/bold]")
        imported = import_clinvar_variants(session, filepath)
        console.print(f"Imported {imported:,} variants")
    else:
        # Check if data file exists
        filepath = config.CLINVAR_DATA_DIR / "variant_summary.txt.gz"
        if filepath.exists():
            console.print(f"[bold]Importing ClinVar variants from existing file...[/bold]")
            imported = import_clinvar_variants(session, filepath)
            console.print(f"Imported {imported:,} variants")
        else:
            console.print("[yellow]No ClinVar data file found. Run with --download first.[/yellow]")

    if do_match:
        console.print(f"\n[bold]Matching pegRNA entries to ClinVar variants...[/bold]")
        matches = match_pegrna_to_clinvar(session, pathogenic_only=pathogenic_only)
        console.print(f"Created {matches:,} matches")

    # Summary
    from database.models import ClinVarVariant, PegRNAClinVarMatch
    from sqlalchemy import func, distinct

    cv_count = session.query(ClinVarVariant).count()
    match_count = session.query(PegRNAClinVarMatch).count()
    matched_entries = session.query(func.count(distinct(PegRNAClinVarMatch.pegrna_entry_id))).scalar()

    console.print(f"\n[bold green]ClinVar Summary:[/bold green]")
    console.print(f"  Variants in DB: {cv_count:,}")
    console.print(f"  Total matches: {match_count:,}")
    console.print(f"  pegRNA entries with matches: {matched_entries:,}")

    session.close()


@app.command(name="fix-sequences")
def fix_sequences(
    paper_id: Optional[int] = typer.Argument(None, help="Fix a specific paper ID (or all if omitted)"),
    limit: int = typer.Option(50, "--limit", "-l", help="Max papers to process"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Show what would be updated without committing"),
):
    """Re-extract and fill missing sequences for papers already in the database.

    Uses expanded column patterns and full-sequence decomposition to recover
    sequences that were missed during initial extraction.
    """
    from database.models import Paper, PegRNAEntry
    from retrieval.supplementary import (
        list_supplementary_files, download_supplementary_file,
        parse_supplementary_tables, find_pegrna_tables,
    )
    from extraction.rule_based import extract_from_multiple_tables
    from extraction.schemas import PegRNAExtracted
    from sqlalchemy import func, case

    session = get_session()

    # Find papers with missing sequences
    if paper_id:
        papers = [(paper_id,)]
    else:
        papers = session.execute(
            session.query(
                Paper.id,
            ).join(PegRNAEntry).group_by(Paper.id).having(
                func.sum(case(
                    (PegRNAEntry.spacer_sequence.is_(None), 1), else_=0
                )) + func.sum(case(
                    (PegRNAEntry.pbs_sequence.is_(None), 1), else_=0
                )) + func.sum(case(
                    (PegRNAEntry.rtt_sequence.is_(None), 1), else_=0
                )) > 0
            ).order_by(
                (func.sum(case(
                    (PegRNAEntry.spacer_sequence.is_(None), 1), else_=0
                )) + func.sum(case(
                    (PegRNAEntry.pbs_sequence.is_(None), 1), else_=0
                )) + func.sum(case(
                    (PegRNAEntry.rtt_sequence.is_(None), 1), else_=0
                ))).desc()
            ).limit(limit)
        ).all()

    if not papers:
        console.print("[green]No papers with missing sequences found[/green]")
        session.close()
        return

    console.print(f"Processing {len(papers)} papers with missing sequences...")

    total_updated = 0
    total_decomposed = 0

    for (pid,) in papers:
        paper = session.query(Paper).get(pid)
        if not paper:
            continue

        entry_count = session.query(PegRNAEntry).filter_by(paper_id=pid).count()
        missing = session.query(
            func.sum(case((PegRNAEntry.spacer_sequence.is_(None), 1), else_=0)),
            func.sum(case((PegRNAEntry.pbs_sequence.is_(None), 1), else_=0)),
            func.sum(case((PegRNAEntry.rtt_sequence.is_(None), 1), else_=0)),
        ).filter_by(paper_id=pid).first()

        console.print(
            f"\n[bold]Paper {pid} (PMID {paper.pmid})[/bold]: "
            f"{entry_count} entries, missing spacer={missing[0]}, pbs={missing[1]}, rtt={missing[2]}"
        )

        # Step A: Re-extract from supplementary files with expanded patterns
        new_entries = []
        if paper.pmcid:
            save_dir = config.RAW_PAPERS_DIR / paper.pmcid
            cached_files = []
            if save_dir.exists():
                for f in save_dir.iterdir():
                    if f.suffix.lower() in ('.xlsx', '.xls', '.csv', '.tsv'):
                        cached_files.append(f)

            if not cached_files:
                supp_files = list_supplementary_files(paper.pmcid)
                tabular = [f for f in supp_files if f["type"] in ("excel", "csv", "tsv")]
                for supp in tabular:
                    filepath = download_supplementary_file(paper.pmcid, supp["filename"])
                    if filepath:
                        cached_files.append(filepath)

            for filepath in cached_files:
                dfs = parse_supplementary_tables(filepath)
                pegrna_tables = find_pegrna_tables(dfs)
                if pegrna_tables:
                    entries = extract_from_multiple_tables(pegrna_tables)
                    new_entries.extend(entries)

        # Match and update existing entries
        if new_entries:
            entry_dicts = [e.to_db_dict() for e in new_entries]
            updated = bulk_update_sequences_for_paper(session, pid, entry_dicts)
            total_updated += updated
            console.print(f"  Updated {updated} entries from re-extraction")

        # Step B: Decompose full_sequence for entries that have it but missing components
        entries_with_full = (
            session.query(PegRNAEntry)
            .filter_by(paper_id=pid)
            .filter(PegRNAEntry.full_sequence.isnot(None))
            .filter(
                (PegRNAEntry.spacer_sequence.is_(None))
                | (PegRNAEntry.pbs_sequence.is_(None))
                | (PegRNAEntry.rtt_sequence.is_(None))
            )
            .all()
        )
        for entry in entries_with_full:
            # Re-validate through schema to trigger decomposition
            try:
                schema = PegRNAExtracted(
                    spacer_sequence=entry.spacer_sequence,
                    pbs_sequence=entry.pbs_sequence,
                    pbs_length=entry.pbs_length,
                    rtt_sequence=entry.rtt_sequence,
                    rtt_length=entry.rtt_length,
                    full_sequence=entry.full_sequence,
                    extension_sequence=None,
                    entry_name=entry.entry_name,
                )
                changed = False
                if schema.spacer_sequence and not entry.spacer_sequence:
                    entry.spacer_sequence = schema.spacer_sequence
                    changed = True
                if schema.pbs_sequence and not entry.pbs_sequence:
                    entry.pbs_sequence = schema.pbs_sequence
                    entry.pbs_length = schema.pbs_length
                    changed = True
                if schema.rtt_sequence and not entry.rtt_sequence:
                    entry.rtt_sequence = schema.rtt_sequence
                    entry.rtt_length = schema.rtt_length
                    changed = True
                if changed:
                    total_decomposed += 1
            except Exception:
                pass

        if not dry_run:
            session.commit()
        else:
            session.rollback()

    console.print(f"\n[bold green]Fix-sequences complete![/bold green]")
    console.print(f"  Updated from re-extraction: {total_updated}")
    console.print(f"  Decomposed from full_sequence: {total_decomposed}")
    if dry_run:
        console.print("[yellow]Dry run — no changes committed[/yellow]")

    session.close()


@app.command(name="annotate-regions")
def annotate_regions(
    organism: Optional[str] = typer.Option(None, help="Filter by organism (e.g., 'Homo sapiens')"),
    gene: Optional[str] = typer.Option(None, help="Annotate a specific gene only"),
    limit: int = typer.Option(0, "--limit", "-l", help="Max genes to annotate (0=all)"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Preview without saving"),
    batch_size: int = typer.Option(100, "--batch-size", help="Genes per Ensembl API batch"),
):
    """Annotate pegRNA entries with exon/intron target region using Ensembl gene structures."""
    from database.models import PegRNAEntry, GeneStructure
    from database.ensembl import batch_fetch_gene_structures, ORGANISM_MAP
    from database.annotate_region import annotate_entries_for_gene
    from sqlalchemy import func, distinct

    session = get_session()

    # Find distinct (gene, organism) pairs needing annotation
    query = (
        session.query(
            PegRNAEntry.target_gene,
            PegRNAEntry.target_organism,
            func.count(PegRNAEntry.id).label("cnt"),
        )
        .filter(
            PegRNAEntry.target_region.is_(None),
            PegRNAEntry.spacer_sequence.isnot(None),
            PegRNAEntry.target_gene.isnot(None),
            PegRNAEntry.target_organism.isnot(None),
        )
        .group_by(PegRNAEntry.target_gene, PegRNAEntry.target_organism)
    )

    if organism:
        query = query.filter(PegRNAEntry.target_organism.ilike(f"%{organism}%"))
    if gene:
        query = query.filter(PegRNAEntry.target_gene == gene)

    gene_org_pairs = query.all()
    console.print(f"Found {len(gene_org_pairs)} gene-organism pairs to annotate")

    if not gene_org_pairs:
        console.print("[green]All entries already annotated[/green]")
        session.close()
        return

    # Group by organism for efficient batch fetching
    org_genes = {}
    for g, o, cnt in gene_org_pairs:
        org_genes.setdefault(o, []).append((g, cnt))

    total_annotated = 0
    total_not_found = 0
    total_genes_processed = 0
    genes_limit = limit if limit > 0 else 999999999

    for org, genes_list in org_genes.items():
        # Check if organism is supported
        org_supported = any(k.lower() == org.lower() for k in ORGANISM_MAP)
        if not org_supported:
            console.print(f"[yellow]Skipping unsupported organism: {org} ({len(genes_list)} genes)[/yellow]")
            continue

        gene_names = [g for g, _ in genes_list]
        total_entries = sum(cnt for _, cnt in genes_list)
        console.print(f"\n[bold]{org}[/bold]: {len(gene_names)} genes, {total_entries:,} entries")

        # Process in batches
        for i in range(0, len(gene_names), batch_size):
            if total_genes_processed >= genes_limit:
                break

            batch = gene_names[i:i + batch_size]
            remaining = int(genes_limit - total_genes_processed)
            if remaining < len(batch):
                batch = batch[:remaining]

            console.print(f"  Fetching gene structures {i+1}-{i+len(batch)} of {len(gene_names)}...")
            structures = batch_fetch_gene_structures(session, batch, org)
            session.commit()

            found = len(structures)
            console.print(f"  Found {found}/{len(batch)} genes in Ensembl")

            for sym in batch:
                gs = structures.get(sym)
                if not gs:
                    continue

                annotated, nf, _ = annotate_entries_for_gene(session, sym, gs, dry_run=dry_run)
                total_annotated += annotated
                total_not_found += nf
                total_genes_processed += 1

            if not dry_run:
                session.commit()
            else:
                session.rollback()

            # Progress
            console.print(
                f"  [green]Annotated: {total_annotated:,}[/green] | "
                f"[yellow]Spacer not found: {total_not_found:,}[/yellow] | "
                f"Genes: {total_genes_processed}"
            )

    console.print(f"\n[bold green]Annotation complete![/bold green]")
    console.print(f"  Total annotated: {total_annotated:,}")
    console.print(f"  Spacer not mapped: {total_not_found:,}")
    console.print(f"  Genes processed: {total_genes_processed}")
    if dry_run:
        console.print("[yellow]Dry run — no changes committed[/yellow]")

    session.close()


@app.command(name="annotate-lof")
def annotate_lof(
    dry_run: bool = typer.Option(False, "--dry-run", help="Preview without saving"),
    force: bool = typer.Option(False, "--force", help="Re-classify already-annotated entries"),
    batch_size: int = typer.Option(5000, "--batch-size", help="Commit batch size"),
    computational: bool = typer.Option(False, "--computational", help="Run Phase 2: RTT-based stop codon detection"),
    fetch_cds: bool = typer.Option(False, "--fetch-cds", help="Fetch CDS sequences from Ensembl before computational analysis"),
    gene: Optional[str] = typer.Option(None, "--gene", help="Target a specific gene (computational mode)"),
    organism: Optional[str] = typer.Option(None, "--organism", help="Filter by organism (computational mode)"),
):
    """Classify pegRNA entries by functional effect (Loss-of-Function).

    Phase 1 (default): HGVS notation (Ter/fs patterns), target region,
    and edit description keywords.

    Phase 2 (--computational): Translate RTT sequences against gene CDS
    to detect premature stop codons and frameshifts in coding exons.
    """
    from database.annotate_lof import annotate_lof_entries, annotate_lof_computational

    Session = init_db(str(config.DATABASE_PATH))
    session = Session()

    # Phase 1: Regex-based classification
    console.print("[bold]Phase 1: Regex-based LoF classification...[/bold]")

    counts = annotate_lof_entries(
        session,
        batch_size=batch_size,
        dry_run=dry_run,
        force=force,
    )

    if not dry_run:
        session.commit()

    classified = counts["total"] - counts["skipped"]
    console.print(f"\n[bold]Phase 1 results:[/bold]")
    console.print(f"  Total entries processed: {counts['total']:,}")
    console.print(f"  Nonsense (stop codon):   {counts['Nonsense']:,}")
    console.print(f"  Frameshift:              {counts['Frameshift']:,}")
    console.print(f"  Splice disruption:       {counts['Splice disruption']:,}")
    console.print(f"  Knockout (annotated):    {counts['Knockout']:,}")
    console.print(f"  [bold]Total classified as LoF: {classified:,}[/bold]")
    console.print(f"  Not LoF (skipped):       {counts['skipped']:,}")

    # Phase 2: Computational RTT-based detection
    if computational:
        # Optionally fetch CDS sequences first
        if fetch_cds:
            from database.ensembl import batch_fetch_cds_sequences
            console.print(f"\n[bold]Fetching CDS sequences from Ensembl...[/bold]")
            fetched = batch_fetch_cds_sequences(session, organism=organism)
            if not dry_run:
                session.commit()
            console.print(f"  Fetched {fetched} CDS sequences")

        console.print(f"\n[bold]Phase 2: Computational RTT-based stop codon detection...[/bold]")

        comp_counts = annotate_lof_computational(
            session,
            batch_size=batch_size,
            dry_run=dry_run,
            gene=gene,
            organism=organism,
        )

        if not dry_run:
            session.commit()

        comp_classified = comp_counts["Nonsense"] + comp_counts["Frameshift"]
        console.print(f"\n[bold]Phase 2 results:[/bold]")
        console.print(f"  Candidates (exon + RTT): {comp_counts['total']:,}")
        console.print(f"  Nonsense (RTT-predicted):{comp_counts['Nonsense']:,}")
        console.print(f"  Frameshift (RTT-pred.):  {comp_counts['Frameshift']:,}")
        console.print(f"  [bold]Total new LoF found:   {comp_classified:,}[/bold]")
        console.print(f"  No CDS available:        {comp_counts['no_cds']:,}")
        console.print(f"  Errors:                  {comp_counts['errors']:,}")
        console.print(f"  No effect detected:      {comp_counts['skipped']:,}")

    if dry_run:
        console.print("\n[yellow]Dry run — no changes committed[/yellow]")

    session.close()


if __name__ == "__main__":
    app()
