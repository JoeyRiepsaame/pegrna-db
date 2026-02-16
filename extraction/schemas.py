"""Pydantic schemas for pegRNA data validation."""
import math
import re
from typing import Optional
from pydantic import BaseModel, field_validator, model_validator


DNA_PATTERN = re.compile(r"^[ACGTUacgtu]+$")


def _coerce_to_str(v):
    """Convert any value to string, handling numpy/pandas types and NaN."""
    if v is None:
        return None
    # Handle NaN (float or numpy)
    try:
        if isinstance(v, float) and math.isnan(v):
            return None
    except (TypeError, ValueError):
        pass
    # Convert numpy scalars to Python natives
    if hasattr(v, 'item'):
        v = v.item()
    if isinstance(v, (int, float)):
        # Convert numeric to string
        if isinstance(v, float) and v == int(v):
            return str(int(v))
        return str(v)
    return str(v) if v is not None else None


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
    extension_sequence: Optional[str] = None  # Combined RTT+PBS from 3' extension column
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

    @field_validator("three_prime_extension", "target_locus", "target_organism",
                     "edit_description", "intended_mutation", "entry_name",
                     "delivery_method", "raw_source_text",
                     mode="before")
    @classmethod
    def coerce_string_fields(cls, v):
        return _coerce_to_str(v)

    @field_validator("spacer_sequence", "pbs_sequence", "rtt_sequence",
                     "linker_sequence", "full_sequence", "nicking_sgrna_seq",
                     "extension_sequence", mode="before")
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
            "replacement": "substitution", "snp": "substitution",
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

    @field_validator("target_organism", mode="before")
    @classmethod
    def normalize_organism(cls, v):
        if v is None:
            return v
        v_str = _coerce_to_str(v)
        if v_str is None:
            return None
        v_clean = v_str.strip().lower()
        bad_values = {"epegrna", "ngrna", "pegrna", "sgrna", "nicking",
                      "none", "n/a", "na", "-", ""}
        if v_clean in bad_values:
            return None
        organism_map = {
            "human": "Homo sapiens", "humans": "Homo sapiens",
            "mouse": "Mus musculus", "mice": "Mus musculus",
            "rat": "Rattus norvegicus", "zebrafish": "Danio rerio",
            "rice": "Oryza sativa", "wheat": "Triticum aestivum",
            "maize": "Zea mays", "corn": "Zea mays",
            "drosophila": "Drosophila melanogaster",
            "pig": "Sus scrofa", "rabbit": "Oryctolagus cuniculus",
        }
        return organism_map.get(v_clean, v_str.strip())

    @model_validator(mode="after")
    def decompose_full_sequence(self):
        """Decompose full_sequence into spacer + extension using scaffold matching."""
        if not self.full_sequence:
            return self
        if self.spacer_sequence and self.pbs_sequence and self.rtt_sequence:
            return self

        # Common SpCas9 scaffold sequences (longest first for best match)
        SCAFFOLDS = [
            "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
            "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
            "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTG",
        ]

        full = self.full_sequence
        for scaffold in SCAFFOLDS:
            idx = full.find(scaffold)
            if idx >= 0:
                # Spacer is before the scaffold (typically 20nt)
                if not self.spacer_sequence and 15 <= idx <= 25:
                    self.spacer_sequence = full[:idx]
                # Extension (RTT+PBS) is after the scaffold
                ext_start = idx + len(scaffold)
                extension = full[ext_start:]
                if extension and not self.extension_sequence:
                    self.extension_sequence = extension
                break
        return self

    @model_validator(mode="after")
    def split_extension_sequence(self):
        """Split extension_sequence (RTT+PBS) into individual components."""
        if not self.extension_sequence:
            return self
        if self.pbs_sequence and self.rtt_sequence:
            return self

        ext = self.extension_sequence
        # If both lengths known, split precisely
        if self.pbs_length and self.rtt_length:
            if self.pbs_length + self.rtt_length == len(ext):
                if not self.rtt_sequence:
                    self.rtt_sequence = ext[:self.rtt_length]
                if not self.pbs_sequence:
                    self.pbs_sequence = ext[self.rtt_length:]
        # If only PBS length known: PBS is at the 3' end
        elif self.pbs_length and not self.pbs_sequence and self.pbs_length < len(ext):
            self.pbs_sequence = ext[-self.pbs_length:]
            if not self.rtt_sequence:
                self.rtt_sequence = ext[:-self.pbs_length]
        # If only RTT length known: RTT is at the 5' end
        elif self.rtt_length and not self.rtt_sequence and self.rtt_length < len(ext):
            self.rtt_sequence = ext[:self.rtt_length]
            if not self.pbs_sequence:
                self.pbs_sequence = ext[self.rtt_length:]
        return self

    @model_validator(mode="after")
    def infer_lengths(self):
        """Infer PBS/RTT lengths from sequences if not provided."""
        if self.pbs_sequence and not self.pbs_length:
            self.pbs_length = len(self.pbs_sequence)
        if self.rtt_sequence and not self.rtt_length:
            self.rtt_length = len(self.rtt_sequence)
        return self

    @model_validator(mode="after")
    def build_extension_from_components(self):
        """Build extension_sequence from RTT + PBS if not already set."""
        if not self.extension_sequence and self.rtt_sequence and self.pbs_sequence:
            self.extension_sequence = self.rtt_sequence + self.pbs_sequence
        return self

    def to_db_dict(self) -> dict:
        """Convert to dict for database insertion."""
        return self.model_dump(exclude_none=False)
