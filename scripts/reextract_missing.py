"""Re-extract missing sequences for papers with incomplete data.

This script:
1. Finds papers with missing spacer/PBS/RTT sequences
2. Re-downloads and re-parses supplementary files with expanded patterns
3. Matches new extractions to existing entries (never creates duplicates)
4. Decomposes full_sequence into components using scaffold matching
5. Reports before/after statistics
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import sqlite3
from datetime import datetime

from rich.console import Console
from rich.table import Table
from sqlalchemy import func, case, distinct

import config
from database.models import init_db, Paper, PegRNAEntry
from database.operations import bulk_update_sequences_for_paper
from retrieval.supplementary import (
    list_supplementary_files, download_supplementary_file,
    parse_supplementary_tables, find_pegrna_tables,
)
from extraction.rule_based import extract_from_multiple_tables
from extraction.schemas import PegRNAExtracted

console = Console()


def get_missing_stats(session):
    """Get current missing sequence statistics."""
    total = session.query(PegRNAEntry).count()
    missing_spacer = session.query(PegRNAEntry).filter(
        PegRNAEntry.spacer_sequence.is_(None)
    ).count()
    missing_pbs = session.query(PegRNAEntry).filter(
        PegRNAEntry.pbs_sequence.is_(None)
    ).count()
    missing_rtt = session.query(PegRNAEntry).filter(
        PegRNAEntry.rtt_sequence.is_(None)
    ).count()
    has_full = session.query(PegRNAEntry).filter(
        PegRNAEntry.full_sequence.isnot(None)
    ).count()
    has_pbs_len_no_seq = session.query(PegRNAEntry).filter(
        PegRNAEntry.pbs_length.isnot(None),
        PegRNAEntry.pbs_sequence.is_(None),
    ).count()

    return {
        "total": total,
        "missing_spacer": missing_spacer,
        "missing_pbs": missing_pbs,
        "missing_rtt": missing_rtt,
        "has_full": has_full,
        "has_pbs_len_no_seq": has_pbs_len_no_seq,
    }


def print_stats(stats, label=""):
    """Print statistics table."""
    table = Table(title=f"Sequence Completeness {label}")
    table.add_column("Metric", style="bold")
    table.add_column("Count", justify="right")
    table.add_column("% of Total", justify="right")

    total = stats["total"]
    for key, val in stats.items():
        if key == "total":
            table.add_row("Total entries", f"{val:,}", "100%")
        else:
            pct = f"{val/total*100:.1f}%" if total > 0 else "0%"
            table.add_row(key.replace("_", " ").title(), f"{val:,}", pct)

    console.print(table)


def get_papers_needing_fix(session, limit=50):
    """Find papers with the most missing sequences."""
    results = (
        session.query(
            Paper.id,
            Paper.pmid,
            Paper.pmcid,
            Paper.title,
            func.count(PegRNAEntry.id).label("total"),
            func.sum(case((PegRNAEntry.spacer_sequence.is_(None), 1), else_=0)).label("miss_spacer"),
            func.sum(case((PegRNAEntry.pbs_sequence.is_(None), 1), else_=0)).label("miss_pbs"),
            func.sum(case((PegRNAEntry.rtt_sequence.is_(None), 1), else_=0)).label("miss_rtt"),
            func.sum(case((PegRNAEntry.full_sequence.isnot(None), 1), else_=0)).label("has_full"),
        )
        .join(PegRNAEntry)
        .group_by(Paper.id)
        .having(
            func.sum(case((PegRNAEntry.spacer_sequence.is_(None), 1), else_=0))
            + func.sum(case((PegRNAEntry.pbs_sequence.is_(None), 1), else_=0))
            + func.sum(case((PegRNAEntry.rtt_sequence.is_(None), 1), else_=0))
            > 0
        )
        .order_by(
            (
                func.sum(case((PegRNAEntry.spacer_sequence.is_(None), 1), else_=0))
                + func.sum(case((PegRNAEntry.pbs_sequence.is_(None), 1), else_=0))
                + func.sum(case((PegRNAEntry.rtt_sequence.is_(None), 1), else_=0))
            ).desc()
        )
        .limit(limit)
        .all()
    )
    return results


def reextract_paper(session, paper, stats_before=None):
    """Re-extract sequences for a single paper."""
    pid = paper.id
    updated_from_supp = 0
    decomposed = 0

    # Step 1: Re-extract from supplementary files
    if paper.pmcid:
        save_dir = config.RAW_PAPERS_DIR / paper.pmcid
        cached_files = []
        if save_dir.exists():
            for f in save_dir.iterdir():
                if f.suffix.lower() in ('.xlsx', '.xls', '.csv', '.tsv'):
                    cached_files.append(f)

        if not cached_files:
            try:
                supp_files = list_supplementary_files(paper.pmcid)
                tabular = [f for f in supp_files if f["type"] in ("excel", "csv", "tsv")]
                for supp in tabular:
                    filepath = download_supplementary_file(paper.pmcid, supp["filename"])
                    if filepath:
                        cached_files.append(filepath)
            except Exception as e:
                console.print(f"  [yellow]Error downloading supps: {e}[/yellow]")

        new_entries = []
        for filepath in cached_files:
            try:
                dfs = parse_supplementary_tables(filepath)
                pegrna_tables = find_pegrna_tables(dfs)
                if pegrna_tables:
                    entries = extract_from_multiple_tables(pegrna_tables)
                    new_entries.extend(entries)
            except Exception as e:
                console.print(f"  [yellow]Error parsing {filepath.name}: {e}[/yellow]")

        if new_entries:
            entry_dicts = [e.to_db_dict() for e in new_entries]
            updated_from_supp = bulk_update_sequences_for_paper(session, pid, entry_dicts)

    # Step 2: Decompose full_sequence entries
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
                decomposed += 1
        except Exception:
            pass

    # Step 3: Process entries with extension_sequence via re-extraction
    # (entries that have pbs_length/rtt_length but no pbs_sequence/rtt_sequence
    #  and an extension column was found during re-extraction)
    entries_with_lengths = (
        session.query(PegRNAEntry)
        .filter_by(paper_id=pid)
        .filter(
            PegRNAEntry.pbs_length.isnot(None),
            PegRNAEntry.rtt_length.isnot(None),
            PegRNAEntry.pbs_sequence.is_(None),
            PegRNAEntry.rtt_sequence.is_(None),
        )
        .all()
    )
    # These entries might be recoverable if their extension was set during
    # the bulk_update step above (via the schema decomposition).
    # The schema model_validator should have already handled this.

    return updated_from_supp, decomposed


def main():
    console.print("[bold]=" * 80)
    console.print("[bold]RE-EXTRACTION: Fixing missing sequences with expanded patterns[/bold]")
    console.print("[bold]=" * 80)

    Session = init_db(str(config.DATABASE_PATH))
    session = Session()

    # Before stats
    console.print("\n[bold]BEFORE:[/bold]")
    before = get_missing_stats(session)
    print_stats(before, "(Before)")

    # Find papers
    papers_data = get_papers_needing_fix(session, limit=150)
    console.print(f"\nFound {len(papers_data)} papers with missing sequences\n")

    total_updated = 0
    total_decomposed = 0

    for i, row in enumerate(papers_data):
        pid, pmid, pmcid, title = row[0], row[1], row[2], row[3]
        total_entries, miss_spacer, miss_pbs, miss_rtt, has_full = (
            row[4], row[5], row[6], row[7], row[8]
        )
        total_missing = miss_spacer + miss_pbs + miss_rtt

        console.print(
            f"[{i+1}/{len(papers_data)}] Paper {pid} (PMID {pmid}): "
            f"{total_entries} entries, {total_missing} missing fields"
        )

        paper = session.query(Paper).get(pid)
        if not paper:
            continue

        updated, decomposed = reextract_paper(session, paper)
        total_updated += updated
        total_decomposed += decomposed

        if updated or decomposed:
            console.print(
                f"  [green]Updated {updated} from supp, "
                f"decomposed {decomposed} from full_seq[/green]"
            )
            session.commit()
        else:
            console.print(f"  [yellow]No updates[/yellow]")

    # After stats
    console.print(f"\n\n[bold]AFTER:[/bold]")
    after = get_missing_stats(session)
    print_stats(after, "(After)")

    # Improvement summary
    console.print(f"\n[bold green]SUMMARY:[/bold green]")
    console.print(f"  Updated from re-extraction: {total_updated}")
    console.print(f"  Decomposed from full_sequence: {total_decomposed}")
    console.print(f"  Spacer: {before['missing_spacer']:,} -> {after['missing_spacer']:,} "
                  f"(recovered {before['missing_spacer'] - after['missing_spacer']:,})")
    console.print(f"  PBS: {before['missing_pbs']:,} -> {after['missing_pbs']:,} "
                  f"(recovered {before['missing_pbs'] - after['missing_pbs']:,})")
    console.print(f"  RTT: {before['missing_rtt']:,} -> {after['missing_rtt']:,} "
                  f"(recovered {before['missing_rtt'] - after['missing_rtt']:,})")

    session.close()


if __name__ == "__main__":
    main()
