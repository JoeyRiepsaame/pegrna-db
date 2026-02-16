"""Database CRUD operations."""
import re
from datetime import datetime
from difflib import SequenceMatcher
from typing import Optional
from sqlalchemy.orm import Session
from sqlalchemy import func
from database.models import Paper, PegRNAEntry


def get_or_create_paper(
    session: Session,
    pmid: Optional[str] = None,
    doi: Optional[str] = None,
    title: Optional[str] = None,
    **kwargs,
) -> tuple[Paper, bool]:
    """Get existing paper or create new one. Returns (paper, created).

    Also checks for title-based duplicates (preprint vs published).
    """
    paper = None
    if pmid:
        paper = session.query(Paper).filter_by(pmid=pmid).first()
    if not paper and doi:
        paper = session.query(Paper).filter_by(doi=doi).first()

    # Title-based dedup for preprint/published pairs
    if not paper and title:
        existing = session.query(Paper).filter(Paper.title.isnot(None)).all()
        for ep in existing:
            sim = SequenceMatcher(
                None, title.lower().strip(), (ep.title or "").lower().strip()
            ).ratio()
            if sim > 0.95:
                paper = ep
                break

    if paper:
        return paper, False

    paper = Paper(pmid=pmid, doi=doi, title=title, **kwargs)
    session.add(paper)
    session.flush()
    return paper, True


def add_pegrna_entry(session: Session, paper_id: int, **kwargs) -> PegRNAEntry:
    """Add a pegRNA entry to the database."""
    entry = PegRNAEntry(paper_id=paper_id, **kwargs)
    session.add(entry)
    session.flush()
    return entry


def bulk_add_entries(session: Session, paper_id: int, entries: list) -> int:
    """Add multiple pegRNA entries. Accepts dicts or PegRNAExtracted objects."""
    count = 0
    for entry_data in entries:
        if hasattr(entry_data, 'to_db_dict'):
            entry_data = entry_data.to_db_dict()
        entry = PegRNAEntry(paper_id=paper_id, **entry_data)
        session.add(entry)
        count += 1
    session.flush()
    return count


def update_paper_status(
    session: Session,
    paper_id: int,
    status: str,
    method: Optional[str] = None,
    notes: Optional[str] = None,
):
    """Update paper extraction status."""
    paper = session.query(Paper).get(paper_id)
    if paper:
        paper.extraction_status = status
        paper.extraction_date = datetime.utcnow()
        if method:
            paper.extraction_method = method
        if notes:
            paper.notes = notes
        session.flush()


def search_entries(
    session: Session,
    target_gene: Optional[str] = None,
    edit_type: Optional[str] = None,
    cell_type: Optional[str] = None,
    prime_editor: Optional[str] = None,
    pegrna_type: Optional[str] = None,
    min_efficiency: Optional[float] = None,
    max_efficiency: Optional[float] = None,
    target_organism: Optional[str] = None,
    editing_technology: Optional[str] = None,
    validated_only: bool = False,
    sort_by: Optional[str] = None,
    sort_desc: bool = False,
    limit: int = 100,
    offset: int = 0,
) -> list[PegRNAEntry]:
    """Search pegRNA entries with filters and sorting."""
    from sqlalchemy import or_
    query = session.query(PegRNAEntry)

    if target_gene:
        query = query.filter(PegRNAEntry.target_gene.ilike(f"%{target_gene}%"))
    if edit_type:
        query = query.filter(PegRNAEntry.edit_type.ilike(f"%{edit_type}%"))
    if cell_type:
        query = query.filter(PegRNAEntry.cell_type.ilike(f"%{cell_type}%"))
    if prime_editor:
        query = query.filter(PegRNAEntry.prime_editor.ilike(f"%{prime_editor}%"))
    if pegrna_type:
        query = query.filter(PegRNAEntry.pegrna_type.ilike(f"%{pegrna_type}%"))
    if min_efficiency is not None:
        query = query.filter(PegRNAEntry.editing_efficiency >= min_efficiency)
    if max_efficiency is not None:
        query = query.filter(PegRNAEntry.editing_efficiency <= max_efficiency)
    if target_organism:
        query = query.filter(PegRNAEntry.target_organism.ilike(f"%{target_organism}%"))
    if editing_technology == "Prime Editing":
        query = query.filter(or_(
            *[PegRNAEntry.prime_editor.ilike(f"%{kw}%") for kw in
              ["PE1", "PE2", "PE3", "PE4", "PE5", "PE6", "PE7", "PEmax", "PECO",
               "twinPE", "mPE", "SaPE", "TJ-PE", "ePE", "pegRNA", "prime"]]
        ))
    elif editing_technology == "Base Editing":
        query = query.filter(or_(
            PegRNAEntry.prime_editor.ilike("%ABE%"),
            PegRNAEntry.prime_editor.ilike("%CBE%"),
        ))
    if validated_only:
        query = query.filter(PegRNAEntry.validated.is_(True))

    # Sorting
    sort_map = {
        "Efficiency": PegRNAEntry.editing_efficiency,
        "Gene": PegRNAEntry.target_gene,
        "Edit Type": PegRNAEntry.edit_type,
        "PE Version": PegRNAEntry.prime_editor,
    }
    if sort_by and sort_by in sort_map:
        col = sort_map[sort_by]
        query = query.order_by(col.desc().nullslast() if sort_desc else col.asc().nullsfirst())
    else:
        query = query.order_by(PegRNAEntry.id)

    return query.offset(offset).limit(limit).all()


def sequence_search(
    session: Session,
    query_sequence: str,
    max_distance: int = 3,
    search_field: str = "spacer",
    limit: int = 100,
) -> list[dict]:
    """Search for pegRNA entries by sequence similarity using edit distance.

    Args:
        query_sequence: DNA sequence to search for (ACGT).
        max_distance: Maximum edit distance (0 = exact, 1-5 = fuzzy).
        search_field: Which sequence field: spacer, pbs, rtt, full.
        limit: Maximum results to return.

    Returns:
        List of dicts with keys: entry, distance.
    """
    import edlib

    query_sequence = query_sequence.strip().upper().replace("U", "T")
    if not re.match(r"^[ACGT]+$", query_sequence):
        return []

    field_map = {
        "spacer": PegRNAEntry.spacer_sequence,
        "pbs": PegRNAEntry.pbs_sequence,
        "rtt": PegRNAEntry.rtt_sequence,
        "full": PegRNAEntry.full_sequence,
    }
    column = field_map.get(search_field, PegRNAEntry.spacer_sequence)

    # Pre-filter by length to reduce comparisons
    qlen = len(query_sequence)
    entries = (
        session.query(PegRNAEntry)
        .filter(column.isnot(None))
        .filter(func.length(column).between(qlen - max_distance, qlen + max_distance))
        .all()
    )

    results = []
    for entry in entries:
        attr = f"{search_field}_sequence" if search_field != "full" else "full_sequence"
        target_seq = getattr(entry, attr)
        if not target_seq:
            continue

        aln = edlib.align(query_sequence, target_seq, mode="NW", task="distance", k=max_distance)
        dist = aln["editDistance"]
        if dist != -1 and dist <= max_distance:
            results.append({"entry": entry, "distance": dist})

    results.sort(key=lambda r: (r["distance"], -(r["entry"].editing_efficiency or 0)))
    return results[:limit]


SEQUENCE_FIELDS = [
    "spacer_sequence", "pbs_sequence", "rtt_sequence",
    "three_prime_extension", "linker_sequence", "full_sequence",
    "nicking_sgrna_seq", "pbs_length", "rtt_length",
    "target_gene", "target_locus", "target_organism",
    "edit_type", "edit_description", "intended_mutation",
    "prime_editor", "cell_type", "delivery_method",
    "editing_efficiency", "product_purity", "indel_frequency",
    "pegrna_type", "entry_name",
]


def _fill_missing_sequences(entry: PegRNAEntry, new_data: dict) -> bool:
    """Fill NULL fields on an existing entry from new_data. Never overwrites.

    Returns True if any field was updated.
    """
    changed = False
    for field in SEQUENCE_FIELDS:
        current = getattr(entry, field, None)
        new_val = new_data.get(field)
        if current is None and new_val is not None:
            setattr(entry, field, new_val)
            changed = True
    return changed


def bulk_update_sequences_for_paper(
    session: Session,
    paper_id: int,
    new_entries: list[dict],
    match_strategy: str = "auto",
) -> int:
    """Match new extraction results to existing entries and fill in NULL sequence fields.

    match_strategy:
    - "name": Match by entry_name
    - "spacer": Match by spacer_sequence
    - "order": Match by row order
    - "auto": Try name, then spacer, then order

    Never creates new rows. Never overwrites existing values.
    """
    existing = (
        session.query(PegRNAEntry)
        .filter_by(paper_id=paper_id)
        .order_by(PegRNAEntry.id)
        .all()
    )

    if not existing:
        return 0

    updated = 0

    if match_strategy in ("auto", "name"):
        name_map = {}
        for entry in existing:
            if entry.entry_name:
                name_map[entry.entry_name.lower().strip()] = entry
        matched_by_name = 0
        for new_data in new_entries:
            name = (new_data.get("entry_name") or "").lower().strip()
            if name and name in name_map:
                if _fill_missing_sequences(name_map[name], new_data):
                    matched_by_name += 1
        if matched_by_name > 0:
            session.flush()
            updated += matched_by_name
            if match_strategy == "auto":
                return updated

    if match_strategy in ("auto", "spacer") and updated == 0:
        spacer_map = {}
        for entry in existing:
            if entry.spacer_sequence:
                spacer_map[entry.spacer_sequence] = entry
        matched_by_spacer = 0
        for new_data in new_entries:
            spacer = (new_data.get("spacer_sequence") or "").strip().upper()
            if spacer and spacer in spacer_map:
                if _fill_missing_sequences(spacer_map[spacer], new_data):
                    matched_by_spacer += 1
        if matched_by_spacer > 0:
            session.flush()
            updated += matched_by_spacer
            if match_strategy == "auto":
                return updated

    if match_strategy in ("auto", "order") and updated == 0:
        matched_by_order = 0
        for i, new_data in enumerate(new_entries):
            if i < len(existing):
                if _fill_missing_sequences(existing[i], new_data):
                    matched_by_order += 1
        if matched_by_order > 0:
            session.flush()
            updated += matched_by_order

    return updated


def get_stats(session: Session) -> dict:
    """Get database statistics."""
    total_papers = session.query(Paper).count()
    completed_papers = session.query(Paper).filter_by(extraction_status="completed").count()
    total_entries = session.query(PegRNAEntry).count()
    validated_entries = session.query(PegRNAEntry).filter_by(validated=True).count()

    # Distinct values
    from sqlalchemy import func, distinct
    genes = session.query(func.count(distinct(PegRNAEntry.target_gene))).scalar()
    cell_types = session.query(func.count(distinct(PegRNAEntry.cell_type))).scalar()
    organisms = session.query(func.count(distinct(PegRNAEntry.target_organism))).scalar()

    return {
        "total_papers": total_papers,
        "completed_papers": completed_papers,
        "total_entries": total_entries,
        "validated_entries": validated_entries,
        "unique_genes": genes,
        "unique_cell_types": cell_types,
        "unique_organisms": organisms,
    }
