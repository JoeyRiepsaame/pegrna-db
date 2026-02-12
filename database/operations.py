"""Database CRUD operations."""
from datetime import datetime
from typing import Optional
from sqlalchemy.orm import Session
from database.models import Paper, PegRNAEntry


def get_or_create_paper(
    session: Session,
    pmid: Optional[str] = None,
    doi: Optional[str] = None,
    **kwargs,
) -> tuple[Paper, bool]:
    """Get existing paper or create new one. Returns (paper, created)."""
    paper = None
    if pmid:
        paper = session.query(Paper).filter_by(pmid=pmid).first()
    if not paper and doi:
        paper = session.query(Paper).filter_by(doi=doi).first()

    if paper:
        return paper, False

    paper = Paper(pmid=pmid, doi=doi, **kwargs)
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
    validated_only: bool = False,
    limit: int = 100,
    offset: int = 0,
) -> list[PegRNAEntry]:
    """Search pegRNA entries with filters."""
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
    if validated_only:
        query = query.filter(PegRNAEntry.validated.is_(True))

    return query.offset(offset).limit(limit).all()


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
