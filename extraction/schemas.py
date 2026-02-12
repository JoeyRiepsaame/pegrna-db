"""Pydantic schemas for pegRNA data validation."""
import re
from typing import Optional
from pydantic import BaseModel, field_validator, model_validator


DNA_PATTERN = re.compile(r"^[ACGTUacgtu]+$")


class PegRNAExtracted(BaseModel):
    """Schema for a single extracted pegRNA entry."""
    model_config = {"coerce_numbers_to_str": True}

    entry_name: Optional[str] = None
    pegrna_type: Optional[str] = None  # pegRNA / epegRNA

    # Sequences
    spacer_sequence: Optional[str] = None
    pbs_sequence: Optional[str] = None
    pbs_length: Optional[int] = None
    rtt_sequence: Optional[str] = None
    rtt_length: Optional[int] = None
    three_prime_extension: Optional[str] = None
    linker_sequence: Optional[str] = None
    full_sequence: Optional[str] = None
    nicking_sgrna_seq: Optional[str] = None

    # Target
    target_gene: Optional[str] = None
    target_locus: Optional[str] = None
    target_organism: Optional[str] = None

    # Edit
    edit_type: Optional[str] = None
    edit_description: Optional[str] = None
    intended_mutation: Optional[str] = None

    # Experimental conditions
    prime_editor: Optional[str] = None
    cell_type: Optional[str] = None
    delivery_method: Optional[str] = None

    # Results
    editing_efficiency: Optional[float] = None
    product_purity: Optional[float] = None
    indel_frequency: Optional[float] = None

    # Extraction metadata
    confidence_score: float = 0.5
    raw_source_text: Optional[str] = None

    @field_validator("spacer_sequence", "pbs_sequence", "rtt_sequence",
                     "linker_sequence", "full_sequence", "nicking_sgrna_seq",
                     mode="before")
    @classmethod
    def clean_sequence(cls, v):
        if v is None:
            return v
        v = str(v).strip().upper().replace(" ", "").replace("-", "")
        # Convert U to T for DNA context
        v = v.replace("U", "T")
        if v and not DNA_PATTERN.match(v):
            return None  # Invalid sequence, discard
        return v if v else None

    @field_validator("editing_efficiency", "product_purity", "indel_frequency",
                     mode="before")
    @classmethod
    def clean_percentage(cls, v):
        if v is None:
            return v
        try:
            val = float(v)
            if val < 0:
                return None
            if val > 100:
                return None
            return val
        except (ValueError, TypeError):
            return None

    @field_validator("pbs_length", "rtt_length", mode="before")
    @classmethod
    def clean_length(cls, v):
        if v is None:
            return v
        try:
            val = int(v)
            return val if 0 < val < 200 else None
        except (ValueError, TypeError):
            return None

    @field_validator("edit_type", mode="before")
    @classmethod
    def normalize_edit_type(cls, v):
        if v is None:
            return v
        v = str(v).lower().strip()
        mapping = {
            "sub": "substitution", "substitution": "substitution",
            "point mutation": "substitution", "snv": "substitution",
            "replacement": "substitution",
            "ins": "insertion", "insertion": "insertion",
            "del": "deletion", "deletion": "deletion",
        }
        return mapping.get(v, v)

    @field_validator("pegrna_type", mode="before")
    @classmethod
    def normalize_pegrna_type(cls, v):
        if v is None:
            return v
        v = str(v).lower().strip()
        if "epegrna" in v or "engineered" in v:
            return "epegRNA"
        if "pegrna" in v or "peg" in v:
            return "pegRNA"
        return v

    @model_validator(mode="after")
    def infer_lengths(self):
        """Infer PBS/RTT lengths from sequences if not provided."""
        if self.pbs_sequence and not self.pbs_length:
            self.pbs_length = len(self.pbs_sequence)
        if self.rtt_sequence and not self.rtt_length:
            self.rtt_length = len(self.rtt_sequence)
        return self

    def to_db_dict(self) -> dict:
        """Convert to dict for database insertion."""
        return self.model_dump(exclude_none=False)
