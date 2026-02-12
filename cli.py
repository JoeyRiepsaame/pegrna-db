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
    search_entries, get_stats,
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
    webapp_path = Path(__file__).parent / "webapp" / "app.py"
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
    from retrieval.pmc_fetcher import fetch_full_text, extract_text_from_bioc, fetch_by_pmid
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
                    for df in dfs[:3]:  # Limit to first 3 sheets
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


if __name__ == "__main__":
    app()
