"""Annotate pegRNA entries with functional effect (Loss-of-Function) classification.

Rule-based classification using HGVS notation in entry_name, target_region,
and edit_description fields. No external API calls needed.

Classification priority:
1. Nonsense (stop codon) - HGVS "Ter" pattern or "nonsense"/"stop codon" keywords
2. Frameshift - HGVS "fs" pattern
3. Splice disruption - target_region == "Splice site"
4. Knockout (annotated) - edit_description contains KO/LoF/knockout keywords
"""
import re
from typing import Optional

from rich.console import Console
from sqlalchemy.orm import Session

from database.models import PegRNAEntry

console = Console()

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
