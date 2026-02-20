"""Annotate pegRNA entries with exon/intron target region information."""
import json
from typing import Optional

from rich.console import Console
from sqlalchemy.orm import Session
from sqlalchemy import func

from database.models import PegRNAEntry, GeneStructure

console = Console()

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


def find_spacer_in_gene(
    spacer: str, gene_sequence: str, gene_start: int
) -> Optional[int]:
    """Find the genomic start coordinate of a spacer within a gene sequence.

    Searches both strands. Returns genomic coordinate or None.
    """
    spacer_upper = spacer.upper()
    gene_upper = gene_sequence.upper()

    # Forward strand
    idx = gene_upper.find(spacer_upper)
    if idx >= 0:
        return gene_start + idx

    # Reverse complement
    rc = reverse_complement(spacer_upper)
    idx = gene_upper.find(rc)
    if idx >= 0:
        return gene_start + idx

    return None


def classify_position(
    match_start: int,
    spacer_length: int,
    exon_ranges: list[tuple[int, int]],
    splice_window: int = 2,
) -> tuple[str, str]:
    """Classify a genomic position as Exon, Intron, or Splice site.

    Returns (category, detail) tuple.
    """
    match_end = match_start + spacer_length - 1
    total_exons = len(exon_ranges)

    if total_exons == 0:
        return "Exon", "Single exon / no annotation"

    in_exons = []
    in_introns = []
    near_splice = False

    # Check exon overlaps
    for i, (ex_start, ex_end) in enumerate(exon_ranges):
        if match_start <= ex_end and match_end >= ex_start:
            in_exons.append(i + 1)

        # Check proximity to exon boundaries (splice sites)
        for boundary in (ex_start, ex_end):
            if (match_start - splice_window <= boundary <= match_end + splice_window):
                # Spacer is near this boundary
                if i > 0 and boundary == ex_start:
                    near_splice = True
                if i < total_exons - 1 and boundary == ex_end:
                    near_splice = True

    # Check intron overlaps (gaps between exons)
    for i in range(total_exons - 1):
        intron_start = exon_ranges[i][1] + 1
        intron_end = exon_ranges[i + 1][0] - 1
        if intron_start <= intron_end and match_start <= intron_end and match_end >= intron_start:
            in_introns.append(i + 1)

    # Classify
    if in_exons and in_introns:
        ex = in_exons[0]
        return "Splice site", f"Exon {ex}/intron junction (of {total_exons} exons)"
    elif near_splice and (in_exons or in_introns):
        if in_exons:
            ex = in_exons[0]
            return "Splice site", f"Near exon {ex} boundary (of {total_exons} exons)"
        else:
            intr = in_introns[0]
            return "Splice site", f"Near intron {intr} boundary (of {total_exons} exons)"
    elif in_exons:
        ex = in_exons[0]
        return "Exon", f"Exon {ex} of {total_exons}"
    elif in_introns:
        intr = in_introns[0]
        return "Intron", f"Intron {intr}-{intr + 1} of {total_exons - 1}"
    else:
        # Outside annotated exons/introns (upstream/downstream/UTR)
        if match_end < exon_ranges[0][0]:
            return "Exon", "5' UTR / upstream"
        elif match_start > exon_ranges[-1][1]:
            return "Exon", "3' UTR / downstream"
        return "Exon", "Flanking region"


def annotate_entries_for_gene(
    session: Session,
    gene_symbol: str,
    gene_structure: GeneStructure,
    dry_run: bool = False,
) -> tuple[int, int, int]:
    """Annotate all entries for a gene. Returns (annotated, not_found, already_done)."""
    exon_ranges = json.loads(gene_structure.exon_coordinates)
    gene_seq = gene_structure.gene_sequence
    gene_start = gene_structure.gene_start

    if not gene_seq or not exon_ranges or gene_start is None:
        return 0, 0, 0

    entries = (
        session.query(PegRNAEntry)
        .filter(
            PegRNAEntry.target_gene == gene_symbol,
            PegRNAEntry.spacer_sequence.isnot(None),
            PegRNAEntry.target_region.is_(None),
        )
        .all()
    )

    annotated = 0
    not_found = 0
    already_done = 0

    for entry in entries:
        spacer = entry.spacer_sequence.strip()
        if len(spacer) < 15:
            not_found += 1
            continue

        pos = find_spacer_in_gene(spacer, gene_seq, gene_start)
        if pos is None:
            not_found += 1
            continue

        category, detail = classify_position(pos, len(spacer), exon_ranges)

        if not dry_run:
            entry.target_region = category
            entry.target_region_detail = detail

        annotated += 1

    return annotated, not_found, already_done
