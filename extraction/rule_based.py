"""Rule-based extraction of pegRNA data from tables and text."""
import re
from typing import Optional

import pandas as pd
from rich.console import Console

import config
from extraction.schemas import PegRNAExtracted

console = Console()


def extract_from_dataframe(df: pd.DataFrame) -> list[PegRNAExtracted]:
    """Extract pegRNA entries from a pandas DataFrame.

    Maps column names to pegRNA fields using fuzzy pattern matching,
    then converts each row to a validated PegRNAExtracted object.
    """
    col_mapping = _map_columns(df)
    if not col_mapping:
        console.print("[yellow]No mappable columns found in table[/yellow]")
        return []

    console.print(f"Column mapping: {col_mapping}")

    # Find the description/name column (first column, often unnamed or has a figure reference)
    desc_col = _find_description_column(df, col_mapping)

    entries = []
    current_section = ""  # Track figure/section headers

    for idx, row in df.iterrows():
        try:
            # Check for section headers (e.g., "Fig. 2a" with no data)
            if desc_col:
                desc_val = row.get(desc_col)
                if pd.notna(desc_val):
                    desc_str = str(desc_val).strip()
                    # Section headers: short, start with "Fig." or similar, no sequence data
                    if (desc_str.startswith("Fig.") or desc_str.startswith("Extended Data"))  \
                            and len(desc_str) < 30:
                        current_section = desc_str
                        continue

            entry_data = {}
            for target_field, source_col in col_mapping.items():
                val = row.get(source_col)
                if pd.notna(val):
                    entry_data[target_field] = val

            # Parse the description column for name, gene, and type info
            if desc_col:
                desc_val = row.get(desc_col)
                if pd.notna(desc_val):
                    desc_str = str(desc_val).strip()
                    parsed = _parse_description(desc_str)
                    # Only fill in fields not already mapped
                    for k, v in parsed.items():
                        if k not in entry_data or not entry_data[k]:
                            entry_data[k] = v

            # Set metadata
            source_file = df.attrs.get("source_file", "")
            source_sheet = df.attrs.get("source_sheet", "")
            raw_text = f"Row {idx} from {source_file}"
            if source_sheet:
                raw_text += f" / {source_sheet}"
            if current_section:
                raw_text += f" ({current_section})"
            entry_data["raw_source_text"] = raw_text
            entry_data["confidence_score"] = 0.8

            entry = PegRNAExtracted(**entry_data)
            if entry.spacer_sequence or entry.target_gene or entry.full_sequence:
                entries.append(entry)

        except Exception as e:
            console.print(f"[yellow]Row {idx} skipped: {e}[/yellow]")

    console.print(f"[green]Extracted {len(entries)} entries from table[/green]")
    return entries


def _find_description_column(df: pd.DataFrame, col_mapping: dict) -> Optional[str]:
    """Find the description/name column (usually first, often unnamed)."""
    mapped_cols = set(col_mapping.values())
    for col in df.columns:
        if col not in mapped_cols:
            # Check if this column has string values (descriptions)
            sample = df[col].dropna().head(10)
            if len(sample) > 0 and sample.apply(lambda x: isinstance(x, str)).mean() > 0.5:
                return col
    return None


def _parse_description(desc: str) -> dict:
    """Parse a pegRNA description string to extract metadata.

    Examples:
    - "DNMT1 +1 FLAG insertion, evopreQ1" -> gene=DNMT1, edit_type=insertion, 3prime=evopreQ1
    - "On-target pegRNA, HEK3 +1 T·A to A·T" -> gene=HEK3, type=pegRNA, edit_type=substitution
    """
    result = {}

    # Entry name is the full description
    result["entry_name"] = desc

    # Detect pegRNA type
    desc_lower = desc.lower()
    if "epegrna" in desc_lower or "engineered" in desc_lower:
        result["pegrna_type"] = "epegRNA"
    elif "pegrna" in desc_lower or "peg rna" in desc_lower:
        result["pegrna_type"] = "pegRNA"
    elif "sgrna" in desc_lower or "sgRNA" in desc:
        # Skip sgRNA-only entries (not pegRNAs)
        pass

    # Detect 3' motif
    motif_patterns = {
        "evopreq1": "evopreQ1", "evopreq": "evopreQ1",
        "tevopreq1": "tevopreQ1", "tevopreq": "tevopreQ1",
        "mpknot": "mpknot",
        "xrrna": "xrRNA",
    }
    for pat, motif in motif_patterns.items():
        if pat in desc_lower:
            result["three_prime_extension"] = motif
            if "pegrna_type" not in result:
                result["pegrna_type"] = "epegRNA"
            break

    # Detect edit type
    if "insertion" in desc_lower or "ins " in desc_lower:
        result["edit_type"] = "insertion"
    elif "deletion" in desc_lower or "del " in desc_lower:
        result["edit_type"] = "deletion"
    elif any(x in desc_lower for x in ["substitution", " to ", "→", ">"]):
        result["edit_type"] = "substitution"

    # Detect common gene names
    known_genes = [
        "HEK3", "HEK4", "FANCF", "EMX1", "RUNX1", "VEGFA", "DNMT1", "RNF2",
        "HBB", "HEXA", "PRNP", "APOE", "PCSK9", "TP53", "BRCA1", "BRCA2",
        "CFTR", "HTT", "POLR2A", "AAVS1", "SEC61B", "PPP1R12C",
    ]
    for gene in known_genes:
        if gene in desc or gene.lower() in desc_lower:
            result["target_gene"] = gene
            break

    return result


def extract_from_text(text: str) -> list[PegRNAExtracted]:
    """Extract pegRNA sequences from free text using regex patterns.

    Looks for common patterns:
    - 5'-SEQUENCE-3' notation
    - pegRNA name followed by sequence
    - Tables embedded in text
    """
    entries = []

    # Pattern: pegRNA-name ... 5'-SEQUENCE-3'
    seq_pattern = re.compile(
        r"(?:pegRNA|epegRNA|peg\s*RNA)[\s\-_]*"
        r"(\S+)?"  # optional name
        r"[:\s]+"
        r"5['\u2019]?[\-\s]*"
        r"([ACGTUacgtu]{15,})"  # sequence (min 15 nt)
        r"[\-\s]*3['\u2019]?",
        re.IGNORECASE,
    )

    for match in seq_pattern.finditer(text):
        name = match.group(1)
        seq = match.group(2).upper().replace("U", "T")
        entry = PegRNAExtracted(
            entry_name=name,
            full_sequence=seq,
            pegrna_type="pegRNA",
            confidence_score=0.6,
            raw_source_text=text[max(0, match.start() - 50) : match.end() + 50],
        )
        entries.append(entry)

    # Pattern: spacer sequence in context of prime editing discussion
    spacer_pattern = re.compile(
        r"spacer[:\s]+['\"]?([ACGTacgt]{18,25})['\"]?",
        re.IGNORECASE,
    )
    pbs_pattern = re.compile(
        r"(?:PBS|primer.binding)[:\s]+['\"]?([ACGTacgt]{8,20})['\"]?",
        re.IGNORECASE,
    )
    rtt_pattern = re.compile(
        r"(?:RTT|RT.template)[:\s]+['\"]?([ACGTacgt]{5,50})['\"]?",
        re.IGNORECASE,
    )

    # Extract editing efficiencies from text
    eff_pattern = re.compile(
        r"editing\s+efficiency\s+(?:of\s+)?(\d+(?:\.\d+)?)\s*%",
        re.IGNORECASE,
    )

    if entries:
        console.print(f"[green]Found {len(entries)} sequences in text[/green]")

    return entries


def _map_columns(df: pd.DataFrame) -> dict[str, str]:
    """Map DataFrame columns to pegRNA entry fields.

    Returns dict mapping target_field -> source_column_name.
    """
    mapping = {}
    cols_lower = {str(c).lower().strip(): c for c in df.columns}

    field_patterns = {
        "spacer_sequence": config.SPACER_PATTERNS,
        "pbs_sequence": config.PBS_PATTERNS,
        "rtt_sequence": config.RTT_PATTERNS,
        "editing_efficiency": config.EFFICIENCY_PATTERNS,
        "edit_type": config.EDIT_TYPE_PATTERNS,
        "cell_type": config.CELL_TYPE_PATTERNS,
        "target_gene": config.GENE_PATTERNS,
        "prime_editor": config.PE_VERSION_PATTERNS,
    }

    # Additional direct mappings
    direct_mappings = {
        "pegrna_type": ["pegrna_type", "peg_type", "grna_type", "type"],
        "entry_name": ["name", "id", "pegrna_name", "pegrna_id", "guide_name"],
        "pbs_length": ["pbs_length", "pbs_len", "pbs length"],
        "rtt_length": ["rtt_length", "rtt_len", "rtt length", "rt_length"],
        "target_organism": ["organism", "species"],
        "target_locus": ["locus", "position", "genomic_position", "coordinates"],
        "edit_description": ["edit_description", "edit", "mutation", "variant"],
        "intended_mutation": ["intended_mutation", "desired_edit", "desired_mutation"],
        "product_purity": ["purity", "product_purity", "correct_edit_pct"],
        "indel_frequency": ["indel", "indel_freq", "indel_frequency", "indels"],
        "delivery_method": ["delivery", "delivery_method", "transfection"],
        "nicking_sgrna_seq": ["nicking", "nick_sgrna", "ngsrna", "pe3_nick"],
        "three_prime_extension": ["3_prime", "3prime", "3' motif", "extension",
                                   "motif", "evopreq", "mpknot", "tevopreq"],
        "linker_sequence": ["linker", "linker_sequence"],
        "full_sequence": ["full_sequence", "full_seq", "pegrna_sequence",
                         "pegrna_seq", "full length"],
    }

    all_patterns = {**field_patterns, **direct_mappings}

    for field, patterns in all_patterns.items():
        if field in mapping:
            continue
        for pat in patterns:
            for col_lower, col_original in cols_lower.items():
                if pat in col_lower and field not in mapping:
                    mapping[field] = col_original
                    break
            if field in mapping:
                break

    return mapping


def extract_from_multiple_tables(
    dfs: list[pd.DataFrame],
) -> list[PegRNAExtracted]:
    """Extract pegRNA data from multiple DataFrames."""
    all_entries = []
    for df in dfs:
        entries = extract_from_dataframe(df)
        all_entries.extend(entries)
    return all_entries
