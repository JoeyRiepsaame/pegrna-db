"""SQLAlchemy ORM models for the pegRNA database."""
from datetime import datetime
from sqlalchemy import (
    Column, Integer, Text, Float, Boolean, DateTime, ForeignKey, Index, create_engine,
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
    __table_args__ = (
        Index("ix_pegrna_target_gene", "target_gene"),
        Index("ix_pegrna_edit_type", "edit_type"),
        Index("ix_pegrna_cell_type", "cell_type"),
        Index("ix_pegrna_prime_editor", "prime_editor"),
        Index("ix_pegrna_target_organism", "target_organism"),
        Index("ix_pegrna_editing_efficiency", "editing_efficiency"),
        Index("ix_pegrna_pegrna_type", "pegrna_type"),
        Index("ix_pegrna_delivery_method", "delivery_method"),
        Index("ix_pegrna_pbs_length", "pbs_length"),
        Index("ix_pegrna_rtt_length", "rtt_length"),
    )

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
    extension_sequence = Column(Text, nullable=True)  # Combined RTT+PBS (3' extension after scaffold)
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


class ClinVarVariant(Base):
    __tablename__ = "clinvar_variants"
    __table_args__ = (
        Index("ix_clinvar_gene", "gene_symbol"),
        Index("ix_clinvar_chr_start", "chromosome", "start"),
        Index("ix_clinvar_significance", "clinical_significance"),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    variation_id = Column(Integer, unique=True, nullable=False)
    allele_id = Column(Integer, nullable=True)
    variant_type = Column(Text, nullable=True)
    name = Column(Text, nullable=True)
    gene_symbol = Column(Text, nullable=True)
    clinical_significance = Column(Text, nullable=True)
    phenotype_list = Column(Text, nullable=True)
    review_status = Column(Text, nullable=True)
    chromosome = Column(Text, nullable=True)
    start = Column(Integer, nullable=True)
    stop = Column(Integer, nullable=True)
    reference_allele = Column(Text, nullable=True)
    alternate_allele = Column(Text, nullable=True)
    assembly = Column(Text, nullable=True)
    hgvs_cdna = Column(Text, nullable=True)
    hgvs_protein = Column(Text, nullable=True)
    last_evaluated = Column(Text, nullable=True)
    rs_dbsnp = Column(Text, nullable=True)

    matches = relationship("PegRNAClinVarMatch", back_populates="clinvar_variant",
                           cascade="all, delete-orphan")


class PegRNAClinVarMatch(Base):
    __tablename__ = "pegrna_clinvar_matches"
    __table_args__ = (
        Index("ix_match_pegrna", "pegrna_entry_id"),
        Index("ix_match_clinvar", "clinvar_variant_id"),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    pegrna_entry_id = Column(Integer, ForeignKey("pegrna_entries.id"), nullable=False)
    clinvar_variant_id = Column(Integer, ForeignKey("clinvar_variants.id"), nullable=False)
    match_type = Column(Text, nullable=False)
    match_confidence = Column(Float, default=0.5)
    notes = Column(Text, nullable=True)

    pegrna_entry = relationship("PegRNAEntry", backref="clinvar_matches")
    clinvar_variant = relationship("ClinVarVariant", back_populates="matches")


def init_db(db_path: str) -> sessionmaker:
    """Initialize database and return a session factory."""
    engine = create_engine(f"sqlite:///{db_path}", echo=False)
    Base.metadata.create_all(engine)

    # Migration: add extension_sequence column if missing
    from sqlalchemy import text, inspect
    insp = inspect(engine)
    if "pegrna_entries" in insp.get_table_names():
        cols = [c["name"] for c in insp.get_columns("pegrna_entries")]
        if "extension_sequence" not in cols:
            with engine.begin() as conn:
                conn.execute(text(
                    "ALTER TABLE pegrna_entries ADD COLUMN extension_sequence TEXT"
                ))
                # Back-fill from existing RTT + PBS
                conn.execute(text("""
                    UPDATE pegrna_entries
                    SET extension_sequence = rtt_sequence || pbs_sequence
                    WHERE rtt_sequence IS NOT NULL
                      AND pbs_sequence IS NOT NULL
                """))

    return sessionmaker(bind=engine)
