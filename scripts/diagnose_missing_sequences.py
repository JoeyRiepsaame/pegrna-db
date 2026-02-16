"""Diagnose missing sequence data by inspecting supplementary files for top papers."""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import sqlite3

import config
from retrieval.supplementary import list_supplementary_files, download_supplementary_file, parse_supplementary_tables

DB_PATH = str(config.DATABASE_PATH)

# All current patterns (for checking what's already mapped)
ALL_MAPPED_PATTERNS = set()
for patterns in [config.SPACER_PATTERNS, config.PBS_PATTERNS, config.RTT_PATTERNS,
                 config.EFFICIENCY_PATTERNS, config.EDIT_TYPE_PATTERNS,
                 config.CELL_TYPE_PATTERNS, config.GENE_PATTERNS, config.PE_VERSION_PATTERNS]:
    ALL_MAPPED_PATTERNS.update(p.lower() for p in patterns)

# Additional patterns from rule_based.py direct_mappings
DIRECT_PATTERN_EXTRAS = [
    "pegrna_type", "peg_type", "grna_type", "type",
    "name", "id", "pegrna_name", "pegrna_id", "guide_name", "finalname",
    "pbslength", "pbs_length", "pbs_len",
    "rtlength", "rtt_length", "rtt_len", "rt_length",
    "organism", "species",
    "locus", "position", "genomic_position", "coordinates",
    "edit_description", "edit", "mutation", "variant", "phenotype",
    "intended_mutation", "desired_edit", "referenceallele", "alternateallele",
    "purity", "product_purity", "correct_edit_pct",
    "indel", "indel_freq", "indel_frequency", "indels", "averageindel",
    "delivery", "delivery_method", "transfection",
    "nicking", "nick_sgrna", "ngsrna", "pe3_nick",
    "3_prime", "3prime", "motif", "evopreq", "mpknot", "tevopreq",
    "linker", "linker_sequence",
    "full_sequence", "full_seq", "pegrna_sequence", "pegrna_seq", "full length",
]
ALL_MAPPED_PATTERNS.update(p.lower() for p in DIRECT_PATTERN_EXTRAS)

DNA_RE = re.compile(r'^[ACGTacgtUu]{8,}$')


def is_dna_column(series: pd.Series, min_ratio: float = 0.3) -> bool:
    """Check if a column predominantly contains DNA sequences."""
    sample = series.dropna().head(50)
    if len(sample) == 0:
        return False
    matches = sample.apply(lambda x: bool(DNA_RE.match(str(x).strip())))
    return matches.mean() >= min_ratio


def is_column_mapped(col_name: str) -> bool:
    """Check if a column name would be matched by current patterns."""
    col_lower = col_name.lower().strip()
    for pat in ALL_MAPPED_PATTERNS:
        if pat in col_lower:
            return True
    return False


def get_papers_with_missing_data(limit: int = 30):
    """Find papers with the most missing sequences."""
    conn = sqlite3.connect(DB_PATH)
    query = """
        SELECT
            p.id, p.pmid, p.pmcid, p.title,
            COUNT(*) as total,
            SUM(CASE WHEN e.spacer_sequence IS NULL THEN 1 ELSE 0 END) as missing_spacer,
            SUM(CASE WHEN e.pbs_sequence IS NULL THEN 1 ELSE 0 END) as missing_pbs,
            SUM(CASE WHEN e.rtt_sequence IS NULL THEN 1 ELSE 0 END) as missing_rtt,
            SUM(CASE WHEN e.full_sequence IS NOT NULL THEN 1 ELSE 0 END) as has_full
        FROM pegrna_entries e
        JOIN papers p ON e.paper_id = p.id
        GROUP BY p.id
        HAVING missing_spacer > 0 OR missing_pbs > 0 OR missing_rtt > 0
        ORDER BY (missing_spacer + missing_pbs + missing_rtt) DESC
        LIMIT ?
    """
    rows = conn.execute(query, (limit,)).fetchall()
    conn.close()
    return rows


def diagnose_paper(paper_id, pmid, pmcid, title):
    """Download and inspect supplementary files for a single paper."""
    result = {
        "paper_id": paper_id,
        "pmid": pmid,
        "pmcid": pmcid,
        "title": title,
        "files": [],
        "unmapped_dna_columns": [],
        "all_columns": {},
    }

    if not pmcid:
        print(f"  No PMCID — skipping supplementary download")
        # Check for cached files
        for pdir in [config.RAW_PAPERS_DIR / (pmid or ""), config.RAW_PAPERS_DIR]:
            if pdir.exists():
                for f in pdir.iterdir():
                    if f.suffix.lower() in ('.xlsx', '.xls', '.csv', '.tsv'):
                        print(f"  Found cached: {f.name}")
        return result

    # Check for already downloaded files
    save_dir = config.RAW_PAPERS_DIR / pmcid
    cached_files = []
    if save_dir.exists():
        for f in save_dir.iterdir():
            if f.suffix.lower() in ('.xlsx', '.xls', '.csv', '.tsv'):
                cached_files.append(f)

    if not cached_files:
        # Download supplementary files
        print(f"  Downloading supplementary files...")
        supp_files = list_supplementary_files(pmcid)
        tabular = [f for f in supp_files if f["type"] in ("excel", "csv", "tsv")]
        if not tabular:
            print(f"  No tabular supplementary files found")
            return result

        for supp in tabular:
            filepath = download_supplementary_file(pmcid, supp["filename"])
            if filepath:
                cached_files.append(filepath)

    # Inspect each file
    for filepath in cached_files:
        file_info = {"filename": filepath.name, "sheets": []}
        print(f"  File: {filepath.name}")

        try:
            if filepath.suffix.lower() in ('.xlsx', '.xls'):
                xl = pd.ExcelFile(filepath)
                for sheet_name in xl.sheet_names:
                    try:
                        df = pd.read_excel(xl, sheet_name=sheet_name, header=None, nrows=100)
                        if df.empty:
                            continue

                        # Try to find the header row
                        from retrieval.supplementary import _find_header_row
                        header_row = _find_header_row(df)
                        df = pd.read_excel(xl, sheet_name=sheet_name, header=header_row, nrows=100)

                        cols = [str(c) for c in df.columns]
                        sheet_info = {
                            "name": sheet_name,
                            "columns": cols,
                            "nrows_sample": len(df),
                        }

                        result["all_columns"][f"{filepath.name}/{sheet_name}"] = cols

                        # Check each column
                        for col in cols:
                            col_str = str(col)
                            if is_dna_column(df[col]):
                                mapped = is_column_mapped(col_str)
                                status = "MAPPED" if mapped else "**UNMAPPED**"
                                avg_len = df[col].dropna().head(20).apply(lambda x: len(str(x).strip())).mean()
                                print(f"    [{sheet_name}] DNA column: '{col_str}' ({status}, avg len ~{avg_len:.0f}nt)")
                                if not mapped:
                                    result["unmapped_dna_columns"].append({
                                        "file": filepath.name,
                                        "sheet": sheet_name,
                                        "column": col_str,
                                        "avg_length": round(avg_len, 1),
                                        "sample": str(df[col].dropna().iloc[0])[:50] if len(df[col].dropna()) > 0 else "",
                                    })
                            elif not is_column_mapped(col_str):
                                # Check if it looks like it could be a header we should know about
                                col_lower = col_str.lower()
                                if any(kw in col_lower for kw in ["seq", "dna", "rna", "oligo", "spacer",
                                                                    "guide", "pbs", "rtt", "template",
                                                                    "primer", "scaffold", "extension"]):
                                    sample_val = str(df[col].dropna().iloc[0])[:80] if len(df[col].dropna()) > 0 else ""
                                    print(f"    [{sheet_name}] Interesting unmapped: '{col_str}' -> {sample_val}")

                        file_info["sheets"].append(sheet_info)
                    except Exception as e:
                        print(f"    [{sheet_name}] Error: {e}")

            elif filepath.suffix.lower() in ('.csv', '.tsv'):
                sep = ',' if filepath.suffix.lower() == '.csv' else '\t'
                df = pd.read_csv(filepath, sep=sep, nrows=100)
                cols = [str(c) for c in df.columns]
                result["all_columns"][filepath.name] = cols

                for col in cols:
                    if is_dna_column(df[col]):
                        mapped = is_column_mapped(str(col))
                        status = "MAPPED" if mapped else "**UNMAPPED**"
                        avg_len = df[col].dropna().head(20).apply(lambda x: len(str(x).strip())).mean()
                        print(f"    DNA column: '{col}' ({status}, avg len ~{avg_len:.0f}nt)")
                        if not mapped:
                            result["unmapped_dna_columns"].append({
                                "file": filepath.name,
                                "sheet": "",
                                "column": str(col),
                                "avg_length": round(avg_len, 1),
                            })

        except Exception as e:
            print(f"  Error reading {filepath.name}: {e}")

        result["files"].append(file_info)

    return result


def main():
    print("=" * 80)
    print("DIAGNOSTIC: Inspecting supplementary files for papers with missing sequences")
    print("=" * 80)

    papers = get_papers_with_missing_data(30)
    print(f"\nFound {len(papers)} papers with missing sequence data\n")

    all_results = []
    all_unmapped = []

    for row in papers:
        paper_id, pmid, pmcid, title = row[0], row[1], row[2], row[3]
        total, miss_spacer, miss_pbs, miss_rtt, has_full = row[4], row[5], row[6], row[7], row[8]
        total_missing = miss_spacer + miss_pbs + miss_rtt

        print(f"\n{'='*60}")
        print(f"Paper ID={paper_id} | PMID={pmid} | PMCID={pmcid}")
        print(f"Title: {(title or '')[:80]}")
        print(f"Entries: {total} | Missing: spacer={miss_spacer}, pbs={miss_pbs}, rtt={miss_rtt} | has_full={has_full}")
        print(f"{'='*60}")

        result = diagnose_paper(paper_id, pmid, pmcid, title)
        result["total_entries"] = total
        result["missing_spacer"] = miss_spacer
        result["missing_pbs"] = miss_pbs
        result["missing_rtt"] = miss_rtt
        result["has_full_sequence"] = has_full
        all_results.append(result)
        all_unmapped.extend(result["unmapped_dna_columns"])

    # Summary
    print(f"\n\n{'='*80}")
    print("SUMMARY: All unmapped DNA columns found across all papers")
    print(f"{'='*80}")
    for item in all_unmapped:
        print(f"  PMID {next((r['pmid'] for r in all_results if any(f['filename'] == item['file'] for f in r['files'])), '?')} | "
              f"{item['file']}/{item['sheet']} | "
              f"Column: '{item['column']}' | "
              f"Avg len: ~{item['avg_length']}nt | "
              f"Sample: {item.get('sample', '')[:40]}")

    # Save report
    report_path = config.DATA_DIR / "diagnosis_report.json"
    with open(report_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nFull report saved to: {report_path}")


if __name__ == "__main__":
    main()
