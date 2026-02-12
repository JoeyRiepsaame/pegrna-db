"""Validation utilities for extracted pegRNA data."""
import re
from rich.console import Console

from extraction.schemas import PegRNAExtracted

console = Console()

DNA_PATTERN = re.compile(r"^[ACGT]+$")

KNOWN_GENES = {
    "HEK3", "HEK4", "FANCF", "EMX1", "RUNX1", "VEGFA", "DNMT1", "RNF2",
    "HBB", "HEXA", "PRNP", "APOE", "PCSK9", "TP53", "BRCA1", "BRCA2",
    "CFTR", "HTT", "SMN1", "SMN2", "GBA", "LRRK2", "SNCA", "APP",
    "POLR2A", "AAVS1",
}

KNOWN_CELL_TYPES = {
    "HEK293T", "HEK293", "HeLa", "K562", "U2OS", "A549",
    "iPSC", "hESC", "Jurkat", "MCF7", "NIH3T3", "N2A", "Neuro2A",
    "primary fibroblasts", "CD34+", "T cells",
}

KNOWN_PE_VERSIONS = {
    "PE1", "PE2", "PE2max", "PE3", "PE3b", "PE3max",
    "PE4", "PE4max", "PE5", "PE5max", "PE6", "PE7",
    "PEmax", "PE2*", "PE-nuclease",
}

KNOWN_3PRIME_MOTIFS = {
    "evopreQ1", "tevopreQ1", "mpknot", "xrRNA", "G-quadruplex",
}


def validate_entry(entry: PegRNAExtracted) -> tuple[bool, list[str]]:
    """Validate a pegRNA entry.

    Returns:
        (is_valid, list_of_warnings)
    """
    warnings = []
    is_valid = True

    # Sequence validation
    for field_name in ("spacer_sequence", "pbs_sequence", "rtt_sequence",
                       "full_sequence", "nicking_sgrna_seq"):
        seq = getattr(entry, field_name)
        if seq and not DNA_PATTERN.match(seq):
            warnings.append(f"{field_name} contains non-ACGT characters: {seq[:20]}...")
            is_valid = False

    # Spacer length (typically 20nt, but can be 17-23)
    if entry.spacer_sequence:
        slen = len(entry.spacer_sequence)
        if slen < 17 or slen > 25:
            warnings.append(f"Unusual spacer length: {slen}nt")

    # PBS length (typically 8-17nt)
    if entry.pbs_length:
        if entry.pbs_length < 5 or entry.pbs_length > 25:
            warnings.append(f"Unusual PBS length: {entry.pbs_length}nt")

    # RTT length (typically 10-50nt)
    if entry.rtt_length:
        if entry.rtt_length < 3 or entry.rtt_length > 100:
            warnings.append(f"Unusual RTT length: {entry.rtt_length}nt")

    # Efficiency checks
    if entry.editing_efficiency is not None:
        if entry.editing_efficiency < 0 or entry.editing_efficiency > 100:
            warnings.append(f"Efficiency out of range: {entry.editing_efficiency}")
            is_valid = False

    # Check if epegRNA has 3' extension
    if entry.pegrna_type == "epegRNA" and not entry.three_prime_extension:
        warnings.append("epegRNA missing 3' extension motif info")

    # Must have at least one meaningful field
    has_data = any([
        entry.spacer_sequence,
        entry.full_sequence,
        entry.target_gene,
        entry.editing_efficiency is not None,
    ])
    if not has_data:
        warnings.append("Entry has no meaningful data (no sequence, gene, or efficiency)")
        is_valid = False

    return is_valid, warnings


def validate_batch(entries: list[PegRNAExtracted]) -> tuple[list[PegRNAExtracted], list[dict]]:
    """Validate a batch of entries.

    Returns:
        (valid_entries, list_of_invalid_with_reasons)
    """
    valid = []
    invalid = []

    for entry in entries:
        is_valid, warnings = validate_entry(entry)
        if is_valid:
            valid.append(entry)
            if warnings:
                console.print(f"[yellow]Warning for {entry.entry_name}: {'; '.join(warnings)}[/yellow]")
        else:
            invalid.append({
                "entry": entry,
                "reasons": warnings,
            })

    console.print(
        f"Validation: [green]{len(valid)} valid[/green], "
        f"[red]{len(invalid)} invalid[/red] out of {len(entries)} total"
    )
    return valid, invalid
