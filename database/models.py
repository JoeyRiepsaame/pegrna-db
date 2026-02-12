"""SQLAlchemy ORM models for the pegRNA database."""
from datetime import datetime
from sqlalchemy import (
    Column, Integer, Text, Float, Boolean, DateTime, ForeignKey, create_engine,
)
from sqlalchemy.orm import declarative_base, relationship, sessionmaker

Base = declarative_base()


class Paper(Base):
    __tablename__ = "papers"

    id = Column(Integer, primary_key=True, autoincrement=True)
    pmid = Column(Text, unique=True, nullable=True)
    pmcid = Column(Text, nullable=True)
    doi = Column(Text, nullable=True)
    title = Column(Text, nullable=True)
    authors = Column(Text, nullable=True)
    journal = Column(Text, nullable=True)
    year = Column(Integer, nullable=True)
    abstract = Column(Text, nullable=True)
    extraction_status = Column(Text, default="pending")  # pending/completed/failed/review
    extraction_date = Column(DateTime, nullable=True)
    extraction_method = Column(Text, nullable=True)  # rule_based/llm/manual
    notes = Column(Text, nullable=True)

    entries = relationship("PegRNAEntry", back_populates="paper", cascade="all, delete-orphan")

    def __repr__(self):
        return f"<Paper(id={self.id}, pmid={self.pmid}, title={self.title!r:.50})>"


class PegRNAEntry(Base):
    __tablename__ = "pegrna_entries"

    id = Column(Integer, primary_key=True, autoincrement=True)
    paper_id = Column(Integer, ForeignKey("papers.id"), nullable=False)

    # Identity
    entry_name = Column(Text, nullable=True)
    pegrna_type = Column(Text, nullable=True)  # pegRNA / epegRNA

    # Sequences
    spacer_sequence = Column(Text, nullable=True)
    pbs_sequence = Column(Text, nullable=True)
    pbs_length = Column(Integer, nullable=True)
    rtt_sequence = Column(Text, nullable=True)
    rtt_length = Column(Integer, nullable=True)
    three_prime_extension = Column(Text, nullable=True)  # evopreQ1, mpknot, tevopreQ1, etc.
    linker_sequence = Column(Text, nullable=True)
    full_sequence = Column(Text, nullable=True)
    nicking_sgrna_seq = Column(Text, nullable=True)

    # Target
    target_gene = Column(Text, nullable=True)
    target_locus = Column(Text, nullable=True)  # genomic coordinates
    target_organism = Column(Text, nullable=True)

    # Edit
    edit_type = Column(Text, nullable=True)  # substitution/insertion/deletion/combination
    edit_description = Column(Text, nullable=True)
    intended_mutation = Column(Text, nullable=True)

    # Experimental conditions
    prime_editor = Column(Text, nullable=True)  # PE2/PE3/PE3b/PE4/PE5/PE7
    cell_type = Column(Text, nullable=True)
    delivery_method = Column(Text, nullable=True)  # plasmid/RNP/mRNA/lentiviral

    # Results
    editing_efficiency = Column(Float, nullable=True)  # percentage
    product_purity = Column(Float, nullable=True)
    indel_frequency = Column(Float, nullable=True)

    # Metadata
    confidence_score = Column(Float, nullable=True)  # 0.0-1.0
    validated = Column(Boolean, default=False)
    raw_source_text = Column(Text, nullable=True)

    paper = relationship("Paper", back_populates="entries")

    def __repr__(self):
        return (
            f"<PegRNAEntry(id={self.id}, gene={self.target_gene}, "
            f"type={self.pegrna_type}, eff={self.editing_efficiency})>"
        )


def init_db(db_path: str) -> sessionmaker:
    """Initialize database and return a session factory."""
    engine = create_engine(f"sqlite:///{db_path}", echo=False)
    Base.metadata.create_all(engine)
    return sessionmaker(bind=engine)
