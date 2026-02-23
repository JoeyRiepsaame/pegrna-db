"""Annotate pegRNA entries with functional effect (Loss-of-Function) classification.

Phase 1 (regex): HGVS notation in entry_name, target_region, edit_description.
Phase 2 (computational): Translate RTT sequences against gene CDS to detect
    premature stop codons and frameshifts.

Classification priority (Phase 1):
1. Nonsense (stop codon) - HGVS "Ter" pattern or "nonsense"/"stop codon" keywords
2. Frameshift - HGVS "fs" pattern
3. Splice disruption - target_region == "Splice site"
4. Knockout (annotated) - edit_description contains KO/LoF/knockout keywords
"""
import json
import re
from typing import Optional

from rich.console import Console
from sqlalchemy.orm import Session

from database.models import PegRNAEntry, GeneStructure

console = Console()

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

# Pre-compiled regex patterns for HGVS protein notation
# Matches p.Xxx123Ter or p.(Xxx123Ter) or p.Xxx123* or p.*123
TER_PATTERN = re.compile(
    r"p\.[\(]?[A-Z][a-z]{2}\d+(?:Ter|\*)[\)]?",
    re.IGNORECASE,
)

# Matches p.Xxx123fs or p.(Xxx123fs) with optional frameTer detail
FS_PATTERN = re.compile(
    r"p\.[\(]?[A-Z][a-z]{2}\d+fs",
    re.IGNORECASE,
)

# Keywords in edit_description that indicate KO/LoF
KO_KEYWORDS = re.compile(
    r"\b(KO|knockout|knock-out|loss[- ]of[- ]function|LoF|gene disruption)\b",
    re.IGNORECASE,
)

NONSENSE_KEYWORDS = re.compile(
    r"\b(nonsense|stop codon|premature stop|premature termination|PTC)\b",
    re.IGNORECASE,
)


def classify_functional_effect(
    entry_name: Optional[str],
    edit_description: Optional[str],
    target_region: Optional[str],
) -> tuple[Optional[str], Optional[str]]:
    """Classify a single entry's functional effect.

    Returns (category, detail) or (None, None) if not LoF.
    Priority: Nonsense > Frameshift > Splice disruption > Knockout.
    """
    combined_text = f"{entry_name or ''} {edit_description or ''}"

    # Rule 1: Nonsense / Stop codon (HGVS Ter pattern)
    ter_match = TER_PATTERN.search(entry_name or "")
    if ter_match:
        return "Nonsense", f"Stop codon: {ter_match.group()}"

    # Also check edit_description for Ter
    ter_match_desc = TER_PATTERN.search(edit_description or "")
    if ter_match_desc:
        return "Nonsense", f"Stop codon: {ter_match_desc.group()}"

    # Nonsense keywords
    nonsense_match = NONSENSE_KEYWORDS.search(combined_text)
    if nonsense_match:
        return "Nonsense", f"Keyword: {nonsense_match.group()}"

    # Rule 2: Frameshift (HGVS fs pattern)
    fs_match = FS_PATTERN.search(entry_name or "")
    if fs_match:
        return "Frameshift", f"Frameshift: {fs_match.group()}"

    fs_match_desc = FS_PATTERN.search(edit_description or "")
    if fs_match_desc:
        return "Frameshift", f"Frameshift: {fs_match_desc.group()}"

    # Rule 3: Splice disruption (from target_region annotation)
    if target_region and target_region.lower() == "splice site":
        return "Splice disruption", "Target region: Splice site"

    # Rule 4: Annotated knockout/LoF
    ko_match = KO_KEYWORDS.search(edit_description or "")
    if ko_match:
        return "Knockout", f"Keyword: {ko_match.group()}"

    ko_match_name = KO_KEYWORDS.search(entry_name or "")
    if ko_match_name:
        return "Knockout", f"Keyword: {ko_match_name.group()}"

    return None, None


def annotate_lof_entries(
    session: Session,
    batch_size: int = 5000,
    dry_run: bool = False,
    force: bool = False,
) -> dict:
    """Annotate all entries that need functional_effect classification.

    Args:
        session: SQLAlchemy session.
        batch_size: Flush every N entries.
        dry_run: Preview without saving.
        force: Re-classify entries that already have a value.

    Returns dict with counts per category plus skipped/total.
    """
    query = session.query(PegRNAEntry)
    if not force:
        query = query.filter(PegRNAEntry.functional_effect.is_(None))

    entries = query.all()

    counts = {
        "Nonsense": 0,
        "Frameshift": 0,
        "Splice disruption": 0,
        "Knockout": 0,
        "skipped": 0,
        "total": len(entries),
    }

    for i, entry in enumerate(entries):
        category, detail = classify_functional_effect(
            entry.entry_name,
            entry.edit_description,
            entry.target_region,
        )

        if category:
            counts[category] += 1
            if not dry_run:
                entry.functional_effect = category
                entry.functional_effect_detail = detail
        else:
            counts["skipped"] += 1

        if not dry_run and (i + 1) % batch_size == 0:
            session.flush()
            console.print(f"  Processed {i + 1:,}/{len(entries):,}...")

    if not dry_run:
        session.flush()

    return counts


# ---------------------------------------------------------------------------
# Phase 2: Computational stop codon / frameshift detection from RTT sequences
# ---------------------------------------------------------------------------

def _reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def _genomic_to_cds_position(
    genomic_pos: int,
    exon_ranges: list[tuple[int, int]],
    strand: int,
) -> Optional[int]:
    """Convert a genomic coordinate to a 0-based CDS position.

    Returns None if the position is not within any exon.
    """
    if strand == 1 or strand is None:
        cds_pos = 0
        for ex_start, ex_end in sorted(exon_ranges):
            if ex_start <= genomic_pos <= ex_end:
                return cds_pos + (genomic_pos - ex_start)
            if genomic_pos > ex_end:
                cds_pos += (ex_end - ex_start + 1)
        return None
    else:  # -1 strand: CDS reads high→low through exons
        cds_pos = 0
        for ex_start, ex_end in sorted(exon_ranges, key=lambda x: x[0], reverse=True):
            if ex_start <= genomic_pos <= ex_end:
                return cds_pos + (ex_end - genomic_pos)
            if genomic_pos < ex_start:
                cds_pos += (ex_end - ex_start + 1)
        return None


def detect_premature_stop(
    rtt_seq: str,
    spacer_seq: str,
    gene_structure: GeneStructure,
) -> tuple[Optional[str], Optional[str]]:
    """Detect if an RTT sequence introduces a premature stop codon or frameshift.

    Algorithm:
    1. Find spacer in genomic sequence → determine strand match
    2. Derive edited coding sequence (RTT or RC(RTT) depending on strand)
    3. Map spacer position to CDS coordinates
    4. Align edited coding seq to CDS window using edlib semi-global
    5. Replace aligned reference with edit, translate, compare proteins

    Returns (category, detail) or (None, None).
    """
    import edlib
    from Bio.Seq import Seq

    cds = gene_structure.cds_sequence
    gene_seq = gene_structure.gene_sequence
    if not cds or not gene_seq:
        return None, None

    cds = cds.upper()
    gene_seq_upper = gene_seq.upper()
    spacer = spacer_seq.strip().upper()
    rtt = rtt_seq.strip().upper().replace("U", "T")

    if len(rtt) < 5 or len(spacer) < 15:
        return None, None

    # Step 1: Find spacer strand
    rc_spacer = _reverse_complement(spacer)
    spacer_on_forward = spacer in gene_seq_upper
    spacer_on_reverse = rc_spacer in gene_seq_upper

    if not spacer_on_forward and not spacer_on_reverse:
        return None, None

    gene_strand = gene_structure.strand or 1
    spacer_strand = 1 if spacer_on_forward else -1

    # Step 2: Derive edited coding sequence
    # RT synthesizes new DNA on the non-target strand = RC(RTT) in 5'→3'.
    # Coding strand edit depends on which strand is coding:
    #   same strand → coding = RC(RTT), different → coding = RTT.
    if spacer_strand == gene_strand:
        edited_coding = _reverse_complement(rtt)
    else:
        edited_coding = rtt

    # Step 3: Map spacer position to CDS
    if spacer_on_forward:
        idx = gene_seq_upper.find(spacer)
    else:
        idx = gene_seq_upper.find(rc_spacer)

    exon_ranges = json.loads(gene_structure.exon_coordinates)

    # Try multiple positions along the spacer (start, mid, end) to find one in an exon
    cds_pos = None
    for offset in [0, len(spacer) // 2, len(spacer) - 1]:
        genomic_pos = gene_structure.gene_start + idx + offset
        cds_pos = _genomic_to_cds_position(genomic_pos, exon_ranges, gene_strand)
        if cds_pos is not None:
            break

    if cds_pos is None or cds_pos < 0 or cds_pos >= len(cds):
        return None, None

    # Step 4: Align edited coding seq to CDS window
    window = 150
    win_start = max(0, cds_pos - window)
    win_end = min(len(cds), cds_pos + window)
    cds_window = cds[win_start:win_end]

    max_dist = max(5, len(edited_coding) // 3)
    aln = edlib.align(edited_coding, cds_window, mode="HW", task="locations", k=max_dist)

    if aln["editDistance"] == -1 or not aln.get("locations"):
        return None, None

    # Prefer alignment where ref length == query length (substitution over indel)
    query_len = len(edited_coding)
    locations = aln["locations"]
    best_loc = min(
        locations,
        key=lambda loc: abs((loc[1] - loc[0] + 1) - query_len),
    )
    win_ref_start = best_loc[0]
    win_ref_end = best_loc[1] + 1  # edlib end is inclusive
    ref_segment = cds_window[win_ref_start:win_ref_end]

    # Map back to full CDS coordinates
    full_ref_start = win_start + win_ref_start
    full_ref_end = win_start + win_ref_end

    # Step 5: Build edited CDS and translate
    edited_cds = cds[:full_ref_start] + edited_coding + cds[full_ref_end:]

    def _translate(seq):
        trimmed = seq[: len(seq) - (len(seq) % 3)]
        if len(trimmed) < 3:
            return ""
        return str(Seq(trimmed).translate())

    ref_protein = _translate(cds)
    edited_protein = _translate(edited_cds)

    if not ref_protein or not edited_protein:
        return None, None

    # Step 6: Compare for premature stop codon
    ref_stop = ref_protein.find("*")
    edited_stop = edited_protein.find("*")

    if edited_stop >= 0 and (ref_stop < 0 or edited_stop < ref_stop):
        aa_pos = edited_stop + 1
        ref_aa_pos = ref_stop + 1 if ref_stop >= 0 else len(ref_protein)
        return "Nonsense", (
            f"RTT-predicted: premature stop at aa {aa_pos} "
            f"(ref stop at aa {ref_aa_pos})"
        )

    # Check for frameshift (length change not divisible by 3)
    len_diff = len(edited_coding) - len(ref_segment)
    if len_diff != 0 and len_diff % 3 != 0:
        direction = "insertion" if len_diff > 0 else "deletion"
        return "Frameshift", (
            f"RTT-predicted: {abs(len_diff)}bp {direction} (frameshift)"
        )

    return None, None


def annotate_lof_computational(
    session: Session,
    batch_size: int = 5000,
    dry_run: bool = False,
    gene: Optional[str] = None,
    organism: Optional[str] = None,
) -> dict:
    """Phase 2: Computational LoF detection from RTT sequences.

    Analyzes entries in coding exons that have RTT + spacer + gene but no
    functional_effect annotation. Never overwrites existing values.

    Auto-fetches missing CDS sequences from Ensembl.
    """
    from database.ensembl import fetch_cds_sequence, ORGANISM_MAP

    # Find candidate entries
    query = (
        session.query(PegRNAEntry)
        .filter(
            PegRNAEntry.functional_effect.is_(None),
            PegRNAEntry.target_region == "Exon",
            PegRNAEntry.rtt_sequence.isnot(None),
            PegRNAEntry.spacer_sequence.isnot(None),
            PegRNAEntry.target_gene.isnot(None),
            PegRNAEntry.target_organism.isnot(None),
        )
    )
    if gene:
        query = query.filter(PegRNAEntry.target_gene == gene)
    if organism:
        query = query.filter(PegRNAEntry.target_organism.ilike(f"%{organism}%"))

    entries = query.all()
    console.print(f"  Candidates: {len(entries):,} entries in coding exons with RTT + spacer")

    # Group by (gene, organism)
    gene_entries: dict[tuple[str, str], list] = {}
    for entry in entries:
        key = (entry.target_gene, entry.target_organism)
        gene_entries.setdefault(key, []).append(entry)

    console.print(f"  Across {len(gene_entries):,} gene-organism pairs")

    counts = {
        "Nonsense": 0,
        "Frameshift": 0,
        "skipped": 0,
        "no_cds": 0,
        "errors": 0,
        "total": len(entries),
    }
    processed_genes = 0

    for (gene_sym, org), group_entries in gene_entries.items():
        # Get gene structure with CDS
        gs = (
            session.query(GeneStructure)
            .filter_by(gene_symbol=gene_sym, organism=org, fetch_status="success")
            .first()
        )
        if not gs:
            counts["no_cds"] += len(group_entries)
            continue

        # Auto-fetch CDS if missing
        if not gs.cds_sequence and gs.ensembl_transcript_id:
            org_info = ORGANISM_MAP.get(org)
            if not org_info:
                for k, v in ORGANISM_MAP.items():
                    if k.lower() == org.lower():
                        org_info = v
                        break
            if org_info:
                _, server = org_info
                cds = fetch_cds_sequence(gs.ensembl_transcript_id, server)
                if cds:
                    gs.cds_sequence = cds
                    session.flush()

        if not gs.cds_sequence:
            counts["no_cds"] += len(group_entries)
            continue

        # Process each entry for this gene
        for entry in group_entries:
            try:
                category, detail = detect_premature_stop(
                    entry.rtt_sequence, entry.spacer_sequence, gs
                )
                if category:
                    counts[category] = counts.get(category, 0) + 1
                    if not dry_run:
                        entry.functional_effect = category
                        entry.functional_effect_detail = detail
                else:
                    counts["skipped"] += 1
            except Exception:
                counts["errors"] += 1

        processed_genes += 1
        if processed_genes % 100 == 0:
            if not dry_run:
                session.flush()
            console.print(
                f"  Processed {processed_genes:,} genes | "
                f"Nonsense: {counts['Nonsense']:,} | "
                f"Frameshift: {counts['Frameshift']:,}"
            )

    if not dry_run:
        session.flush()

    return counts
