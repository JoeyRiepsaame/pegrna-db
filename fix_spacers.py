"""Fix missing spacer sequences for papers with non-standard column formats.

Run: python fix_spacers.py
"""
import ast
import re
import sys
from pathlib import Path

import pandas as pd
from rich.console import Console

sys.path.insert(0, str(Path(__file__).parent))
import config
from database.models import init_db, Paper, PegRNAEntry

console = Console()


def fix_pridict_spacers(session, paper_id: int):
    """Fix PRIDICT Library_1 spacers by reconstructing from protobase columns."""
    f = Path("data/raw_papers/PMC7614945/EMS157167-supplement-Supplementary_Table_6.xlsx")
    if not f.exists():
        console.print(f"[red]PRIDICT file not found: {f}[/red]")
        return 0

    console.print("[cyan]Reading PRIDICT Library_1...[/cyan]")
    df = pd.read_excel(f, sheet_name="Library_1")
    console.print(f"  {len(df)} rows")

    # Find protobase columns
    pb_cols = []
    for i in range(1, 30):
        candidates = [c for c in df.columns if str(c).lower().strip() == f"protobase_{i}"]
        if candidates:
            pb_cols.append(candidates[0])
        else:
            break
    console.print(f"  Found {len(pb_cols)} protobase columns")

    if len(pb_cols) < 15:
        console.print("[red]Not enough protobase columns[/red]")
        return 0

    # Also check wide_initial_target + protospacerlocation_only_initial as fallback
    wide_col = None
    loc_col = None
    for c in df.columns:
        cl = str(c).lower().strip()
        if cl == "wide_initial_target":
            wide_col = c
        elif cl == "protospacerlocation_only_initial":
            loc_col = c

    # Get existing Library_1 entries (ordered by ID to match row order)
    entries = (
        session.query(PegRNAEntry)
        .filter_by(paper_id=paper_id)
        .filter(PegRNAEntry.raw_source_text.like("%Library_1%"))
        .order_by(PegRNAEntry.id)
        .all()
    )
    console.print(f"  {len(entries)} existing Library_1 entries")

    updated = 0
    for i, entry in enumerate(entries):
        if entry.spacer_sequence:
            continue
        if i >= len(df):
            break

        row = df.iloc[i]

        # Method 1: Reconstruct from protobases
        bases = []
        for col in pb_cols:
            val = row.get(col)
            if pd.notna(val):
                b = str(val).strip().upper()
                if len(b) == 1 and b in "ACGT":
                    bases.append(b)
                else:
                    bases = []
                    break
            else:
                bases = []
                break

        spacer = "".join(bases) if bases else None

        # Method 2: Extract from wide_initial_target
        if not spacer and wide_col and loc_col:
            wide_val = row.get(wide_col)
            loc_val = row.get(loc_col)
            if pd.notna(wide_val) and pd.notna(loc_val):
                try:
                    wide_seq = str(wide_val).strip().upper()
                    if isinstance(loc_val, str):
                        loc = ast.literal_eval(loc_val)
                    elif isinstance(loc_val, (list, tuple)):
                        loc = list(loc_val)
                    else:
                        loc = None
                    if isinstance(loc, list) and len(loc) == 2:
                        start, end = int(loc[0]), int(loc[1])
                        if 0 <= start < end <= len(wide_seq):
                            candidate = wide_seq[start:end]
                            if re.match(r"^[ACGT]+$", candidate):
                                spacer = candidate
                except (ValueError, TypeError, SyntaxError):
                    pass

        if spacer and 15 <= len(spacer) <= 25:
            entry.spacer_sequence = spacer
            updated += 1

    if updated > 0:
        session.flush()
    console.print(f"[green]Updated {updated} PRIDICT Library_1 spacers[/green]")
    return updated


def fix_paper_spacers_with_header_offset(
    session, paper_id: int, file_path: str, sheet_name: str,
    header_row: int, source_text_pattern: str,
):
    """Fix spacers for papers where the header row was misread."""
    f = Path(file_path)
    if not f.exists():
        console.print(f"[red]File not found: {f}[/red]")
        return 0

    console.print(f"[cyan]Reading {f.name} / {sheet_name} (header={header_row})...[/cyan]")
    df = pd.read_excel(f, sheet_name=sheet_name, header=header_row)
    console.print(f"  {len(df)} rows, {len(df.columns)} cols")

    # Use updated column mapping
    from extraction.rule_based import _map_columns, _normalize_col_name
    mapping = _map_columns(df)
    console.print(f"  Mapping: {mapping}")

    spacer_col = mapping.get("spacer_sequence")
    if not spacer_col:
        console.print("[yellow]No spacer column found[/yellow]")
        return 0

    # Build spacer lookup from new extraction
    spacers_by_row = {}
    for idx, row in df.iterrows():
        val = row.get(spacer_col)
        if pd.notna(val):
            spacer = str(val).strip().upper().replace(" ", "").replace("-", "")
            # Strip leading lowercase 'g' prefix (sometimes added for transcription start)
            if spacer.startswith("G") and len(spacer) == 21:
                spacer_clean = spacer[1:]
            else:
                spacer_clean = spacer
            if re.match(r"^[ACGT]+$", spacer_clean) and len(spacer_clean) >= 15:
                spacers_by_row[idx] = spacer_clean

    console.print(f"  Found {len(spacers_by_row)} spacers in re-read data")

    # Get existing entries for this sheet
    entries = (
        session.query(PegRNAEntry)
        .filter_by(paper_id=paper_id)
        .filter(PegRNAEntry.raw_source_text.like(f"%{source_text_pattern}%"))
        .filter(PegRNAEntry.spacer_sequence.is_(None))
        .order_by(PegRNAEntry.id)
        .all()
    )
    console.print(f"  {len(entries)} entries missing spacers")

    # Also get entry_name mapping for name-based matching
    entry_name_col = mapping.get("entry_name")

    # Try order-based matching
    updated = 0
    for i, entry in enumerate(entries):
        if i in spacers_by_row:
            entry.spacer_sequence = spacers_by_row[i]
            updated += 1

    if updated > 0:
        session.flush()
    console.print(f"[green]Updated {updated} spacers for {f.name}[/green]")
    return updated


def fix_paper_spacers_reextract(
    session, paper_id: int, file_path: str, sheet_name: str,
    source_text_pattern: str, header_row: int = 0,
):
    """Re-extract spacers using updated column mapping and fill into existing entries."""
    f = Path(file_path)
    if not f.exists():
        console.print(f"[red]File not found: {f}[/red]")
        return 0

    console.print(f"[cyan]Re-extracting {f.name} / {sheet_name}...[/cyan]")
    df = pd.read_excel(f, sheet_name=sheet_name, header=header_row)
    df.attrs["source_file"] = f.name
    df.attrs["source_sheet"] = sheet_name

    from extraction.rule_based import extract_from_dataframe
    new_entries = extract_from_dataframe(df)
    console.print(f"  Re-extracted {len(new_entries)} entries")

    new_with_spacer = sum(1 for e in new_entries if e.spacer_sequence)
    console.print(f"  {new_with_spacer} with spacer")

    if new_with_spacer == 0:
        return 0

    # Get existing entries
    from database.operations import bulk_update_sequences_for_paper
    new_dicts = [e.to_db_dict() for e in new_entries]
    updated = bulk_update_sequences_for_paper(session, paper_id, new_dicts, match_strategy="auto")
    console.print(f"[green]Updated {updated} entries for {f.name}[/green]")
    return updated


def fix_moesm7_by_name_lookup(session, paper_id: int):
    """Fix PMC12062269 MOESM7 entries by matching pegRNA names to MOESM3/MOESM5 sequences."""
    # Build name -> spacer lookup from MOESM3 and MOESM5
    spacer_lookup = {}

    # MOESM3 (pegRNA and ngRNA sequences)
    f3 = Path("data/raw_papers/PMC12062269/41467_2025_59708_MOESM3_ESM.xlsx")
    if f3.exists():
        df3 = pd.read_excel(f3, sheet_name="pegRNA and ngRNA sequences", header=1)
        spacer_col = "Spacer  (5'-3')"
        target_col = "Target"
        if spacer_col in df3.columns and target_col in df3.columns:
            for _, row in df3.iterrows():
                spacer = row.get(spacer_col)
                target = row.get(target_col)
                if pd.notna(spacer) and pd.notna(target):
                    s = str(spacer).strip().upper().replace(" ", "")
                    if s.startswith("G") and len(s) == 21:
                        s = s[1:]
                    if re.match(r"^[ACGT]+$", s) and len(s) >= 15:
                        spacer_lookup[str(target).strip()] = s

    # MOESM5 (CRISPResso2 parameters)
    f5 = Path("data/raw_papers/PMC12062269/41467_2025_59708_MOESM5_ESM.xlsx")
    if f5.exists():
        df5 = pd.read_excel(f5, sheet_name="CRISPResso2 parameters")
        if "prime_editing_pegRNA_spacer_seq" in df5.columns and "target_name" in df5.columns:
            for _, row in df5.iterrows():
                spacer = row.get("prime_editing_pegRNA_spacer_seq")
                target = row.get("target_name")
                if pd.notna(spacer) and pd.notna(target):
                    s = str(spacer).strip().upper()
                    if re.match(r"^[ACGT]+$", s) and len(s) >= 15:
                        spacer_lookup[str(target).strip()] = s

    console.print(f"  Built lookup with {len(spacer_lookup)} target->spacer mappings")

    # Get MOESM7 entries without spacers
    entries = (
        session.query(PegRNAEntry)
        .filter_by(paper_id=paper_id)
        .filter(PegRNAEntry.raw_source_text.like("%MOESM7%"))
        .filter(PegRNAEntry.spacer_sequence.is_(None))
        .all()
    )
    console.print(f"  {len(entries)} MOESM7 entries missing spacers")

    updated = 0
    for entry in entries:
        # Try matching by target_gene
        if entry.target_gene and entry.target_gene in spacer_lookup:
            entry.spacer_sequence = spacer_lookup[entry.target_gene]
            updated += 1

    if updated > 0:
        session.flush()
    console.print(f"[green]Updated {updated} MOESM7 entries by name lookup[/green]")
    return updated


def main():
    Session = init_db(str(config.DATABASE_PATH))
    session = Session()

    total_updated = 0

    try:
        # 1. PRIDICT (PMID 36646933, paper_id=3) - 92K missing
        console.print("\n[bold]=== PRIDICT (PMID 36646933) ===[/bold]")
        paper = session.query(Paper).filter_by(pmid="36646933").first()
        if paper:
            n = fix_pridict_spacers(session, paper.id)
            total_updated += n

        # 2. PMC11754097 (PMID 38987629) - 5.8K missing
        console.print("\n[bold]=== PMID 38987629 ===[/bold]")
        paper = session.query(Paper).filter_by(pmid="38987629").first()
        if paper:
            n = fix_paper_spacers_reextract(
                session, paper.id,
                "data/raw_papers/PMC11754097/41551_2024_1233_MOESM3_ESM.xlsx",
                "ST2 gRNA sequences + SRA files",
                "MOESM3",
            )
            total_updated += n

        # 3. PMC12062269 (PMID 40341582) - 3.5K missing
        console.print("\n[bold]=== PMID 40341582 ===[/bold]")
        paper = session.query(Paper).filter_by(pmid="40341582").first()
        if paper:
            # Try MOESM3 with header=1
            n = fix_paper_spacers_with_header_offset(
                session, paper.id,
                "data/raw_papers/PMC12062269/41467_2025_59708_MOESM3_ESM.xlsx",
                "pegRNA and ngRNA sequences",
                header_row=1,
                source_text_pattern="MOESM3",
            )
            total_updated += n

            # Try MOESM5 re-extract
            n = fix_paper_spacers_reextract(
                session, paper.id,
                "data/raw_papers/PMC12062269/41467_2025_59708_MOESM5_ESM.xlsx",
                "CRISPResso2 parameters",
                "MOESM5",
            )
            total_updated += n

            # Try MOESM7 name lookup
            n = fix_moesm7_by_name_lookup(session, paper.id)
            total_updated += n

        # 4. PMC10869272 (PMID 37142705) - 1.5K missing
        console.print("\n[bold]=== PMID 37142705 ===[/bold]")
        paper = session.query(Paper).filter_by(pmid="37142705").first()
        if paper:
            n = fix_paper_spacers_with_header_offset(
                session, paper.id,
                "data/raw_papers/PMC10869272/41587_2023_1758_MOESM3_ESM.xlsx",
                "SI Table 2 peRNA and AAV design",
                header_row=1,
                source_text_pattern="SI Table 2",
            )
            total_updated += n

        # 5. Other papers - try re-extract with updated column mapping
        console.print("\n[bold]=== Other papers with missing spacers ===[/bold]")
        other_pmids = ["35332138", "38549229", "39737993", "36658146",
                       "41345100", "33927418", "35351879", "38191664",
                       "35436085", "39747083"]
        for pmid in other_pmids:
            paper = session.query(Paper).filter_by(pmid=pmid).first()
            if not paper:
                continue
            missing = (
                session.query(PegRNAEntry)
                .filter_by(paper_id=paper.id)
                .filter(PegRNAEntry.spacer_sequence.is_(None))
                .count()
            )
            if missing == 0:
                continue

            # Find supplementary files
            pmcid = paper.pmcid
            if not pmcid:
                continue
            paper_dir = Path(f"data/raw_papers/{pmcid}")
            if not paper_dir.exists():
                continue

            xlsx_files = list(paper_dir.glob("*.xlsx")) + list(paper_dir.glob("*.xls"))
            for xlsx in xlsx_files:
                try:
                    xl = pd.ExcelFile(xlsx)
                    for sheet in xl.sheet_names:
                        df = pd.read_excel(xl, sheet_name=sheet, nrows=5)
                        from extraction.rule_based import _map_columns
                        mapping = _map_columns(df)
                        if "spacer_sequence" in mapping:
                            console.print(f"  PMID {pmid}: re-extracting {xlsx.name}/{sheet}")
                            n = fix_paper_spacers_reextract(
                                session, paper.id,
                                str(xlsx), sheet,
                                xlsx.name,
                            )
                            total_updated += n
                except Exception as e:
                    console.print(f"  [yellow]Error with {xlsx.name}: {e}[/yellow]")

        # Commit all changes
        session.commit()
        console.print(f"\n[bold green]Total spacers updated: {total_updated}[/bold green]")

        # Print final statistics
        total = session.query(PegRNAEntry).count()
        with_spacer = session.query(PegRNAEntry).filter(
            PegRNAEntry.spacer_sequence.isnot(None)
        ).count()
        console.print(f"Database: {total} entries, {with_spacer} with spacer ({100*with_spacer/total:.1f}%)")

    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        session.rollback()
        raise
    finally:
        session.close()


if __name__ == "__main__":
    main()
