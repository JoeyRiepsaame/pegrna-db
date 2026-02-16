"""ClinVar data download, import, and pegRNA matching."""
import gzip
import csv
import re
from pathlib import Path
from typing import Optional

from rich.console import Console
from sqlalchemy import func, distinct
from sqlalchemy.orm import Session

from database.models import ClinVarVariant, PegRNAClinVarMatch, PegRNAEntry

console = Console()


def download_clinvar(data_dir: Path, url: str) -> Path:
    """Download ClinVar variant_summary.txt.gz from NCBI FTP.

    Returns path to the downloaded file.
    """
    import urllib.request

    data_dir.mkdir(parents=True, exist_ok=True)
    output_path = data_dir / "variant_summary.txt.gz"

    console.print(f"Downloading ClinVar data to {output_path}...")
    urllib.request.urlretrieve(url, str(output_path))
    console.print(f"[green]Downloaded {output_path.stat().st_size / 1e6:.1f} MB[/green]")

    return output_path


def import_clinvar_variants(
    session: Session,
    filepath: Path,
    batch_size: int = 10000,
) -> int:
    """Parse ClinVar variant_summary.txt.gz and import variants matching DB genes.

    Only imports variants for:
    - Genes that exist in the pegRNA database
    - Human assembly (GRCh37 or GRCh38)
    """
    # Get genes in our DB for filtering
    db_genes = set()
    gene_rows = (
        session.query(distinct(PegRNAEntry.target_gene))
        .filter(PegRNAEntry.target_gene.isnot(None))
        .all()
    )
    for row in gene_rows:
        gene = row[0].strip().upper()
        db_genes.add(gene)
    console.print(f"Found {len(db_genes)} unique genes in pegRNA database")

    if not db_genes:
        console.print("[yellow]No genes in database, skipping ClinVar import[/yellow]")
        return 0

    # Get existing variation_ids to avoid duplicates
    existing_ids = set()
    for row in session.query(ClinVarVariant.variation_id).all():
        existing_ids.add(row[0])
    console.print(f"Already have {len(existing_ids)} ClinVar variants in DB")

    imported = 0
    skipped = 0
    batch = []

    console.print("Parsing ClinVar variant_summary.txt.gz...")
    with gzip.open(str(filepath), "rt", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            # Filter: human assembly only
            assembly = row.get("Assembly", "")
            if assembly not in ("GRCh37", "GRCh38"):
                continue

            # Filter: gene must be in our DB
            gene = (row.get("GeneSymbol") or "").strip().upper()
            if gene not in db_genes:
                continue

            # Skip if already imported
            try:
                var_id = int(row.get("VariationID", 0))
            except (ValueError, TypeError):
                continue

            if var_id in existing_ids:
                skipped += 1
                continue

            # Parse fields
            try:
                start = int(row.get("Start", 0)) if row.get("Start") else None
                stop = int(row.get("Stop", 0)) if row.get("Stop") else None
                allele_id = int(row.get("AlleleID", 0)) if row.get("AlleleID") else None
            except (ValueError, TypeError):
                start = stop = allele_id = None

            variant = ClinVarVariant(
                variation_id=var_id,
                allele_id=allele_id,
                variant_type=row.get("Type"),
                name=row.get("Name"),
                gene_symbol=row.get("GeneSymbol"),
                clinical_significance=row.get("ClinicalSignificance") or row.get("ClinSigSimple"),
                phenotype_list=row.get("PhenotypeList"),
                review_status=row.get("ReviewStatus"),
                chromosome=row.get("Chromosome"),
                start=start,
                stop=stop,
                reference_allele=row.get("ReferenceAllele"),
                alternate_allele=row.get("AlternateAllele"),
                assembly=assembly,
                hgvs_cdna=row.get("HGVS(c)"),
                hgvs_protein=row.get("HGVS(p)"),
                last_evaluated=row.get("LastEvaluated"),
                rs_dbsnp=row.get("RS# (dbSNP)"),
            )
            batch.append(variant)
            existing_ids.add(var_id)

            if len(batch) >= batch_size:
                session.bulk_save_objects(batch)
                session.commit()
                imported += len(batch)
                console.print(f"  Imported {imported:,} variants...")
                batch = []

    # Final batch
    if batch:
        session.bulk_save_objects(batch)
        session.commit()
        imported += len(batch)

    console.print(f"[green]Imported {imported:,} ClinVar variants ({skipped:,} already existed)[/green]")
    return imported


def match_pegrna_to_clinvar(
    session: Session,
    pathogenic_only: bool = True,
    max_variants_per_entry: int = 10,
) -> int:
    """Match pegRNA entries to ClinVar variants by gene name.

    For each pegRNA entry, finds the top ClinVar variants (by confidence) for
    the same gene. Only stores matches with confidence > 0.3 (i.e., at least
    some edit-type or mutation alignment beyond just gene name).

    Confidence scoring:
    - 0.3 base (gene name match)
    - +0.2 edit type alignment (substitution->SNV, deletion->Deletion, insertion->Insertion)
    - +0.2 mutation description overlap
    - +0.1 edit description word overlap

    Only matches above 0.3 base are stored, and limited to top N per entry.
    Processes per-gene for memory efficiency.
    """
    # Get genes that have both pegRNA entries and ClinVar variants
    pegrna_genes = set()
    for row in session.query(distinct(PegRNAEntry.target_gene)).filter(PegRNAEntry.target_gene.isnot(None)).all():
        pegrna_genes.add(row[0].strip())

    clinvar_genes = set()
    for row in session.query(distinct(ClinVarVariant.gene_symbol)).filter(ClinVarVariant.gene_symbol.isnot(None)).all():
        clinvar_genes.add(row[0].strip())

    # Case-insensitive matching
    pegrna_gene_map = {g.upper(): g for g in pegrna_genes}
    clinvar_gene_map = {g.upper(): g for g in clinvar_genes}

    common_genes = set(pegrna_gene_map.keys()) & set(clinvar_gene_map.keys())
    console.print(f"Found {len(common_genes)} genes in common between pegRNA DB and ClinVar")

    if not common_genes:
        console.print("[yellow]No matching genes found[/yellow]")
        return 0

    # Clear existing matches
    existing_count = session.query(PegRNAClinVarMatch).count()
    if existing_count > 0:
        console.print(f"Clearing {existing_count:,} existing matches...")
        session.query(PegRNAClinVarMatch).delete()
        session.commit()

    # Edit type mapping for confidence scoring
    EDIT_TYPE_MAP = {
        "substitution": ["single nucleotide variant", "snv", "missense", "nonsense", "synonymous"],
        "deletion": ["deletion", "del", "microdeletion"],
        "insertion": ["insertion", "ins", "duplication", "dup"],
    }

    total_matches = 0

    for i, gene_upper in enumerate(sorted(common_genes)):
        # Get pegRNA entries for this gene
        entries = (
            session.query(PegRNAEntry)
            .filter(func.upper(PegRNAEntry.target_gene) == gene_upper)
            .all()
        )

        # Get ClinVar variants for this gene
        cv_query = session.query(ClinVarVariant).filter(
            func.upper(ClinVarVariant.gene_symbol) == gene_upper
        )
        if pathogenic_only:
            cv_query = cv_query.filter(
                ClinVarVariant.clinical_significance.ilike("%pathogenic%")
            )
        variants = cv_query.all()

        if not entries or not variants:
            continue

        gene_matches = 0
        for entry in entries:
            # Score all variants for this entry, keep top N above threshold
            scored = []
            entry_edit = (entry.edit_type or "").lower()
            entry_mutation = (entry.intended_mutation or "").lower()
            entry_desc = (entry.edit_description or "").lower()

            for variant in variants:
                confidence = 0.3  # Base: gene match

                # Edit type alignment
                variant_type = (variant.variant_type or "").lower()
                variant_name = (variant.name or "").lower()

                for edit_cat, keywords in EDIT_TYPE_MAP.items():
                    if edit_cat in entry_edit:
                        if any(kw in variant_type or kw in variant_name for kw in keywords):
                            confidence += 0.2
                        break

                # Mutation description overlap
                variant_hgvs = ((variant.hgvs_cdna or "") + " " + (variant.hgvs_protein or "")).lower()
                if entry_mutation and variant_hgvs:
                    entry_tokens = set(re.findall(r'[a-z]+\d+[a-z]+|[acgt]>\s*[acgt]|\d+', entry_mutation))
                    variant_tokens = set(re.findall(r'[a-z]+\d+[a-z]+|[acgt]>\s*[acgt]|\d+', variant_hgvs))
                    if entry_tokens and variant_tokens:
                        overlap = len(entry_tokens & variant_tokens) / max(len(entry_tokens), 1)
                        confidence += min(0.2, overlap * 0.2)

                # Edit description word overlap
                phenotype = (variant.phenotype_list or "").lower()
                if entry_desc and phenotype:
                    desc_words = set(entry_desc.split()) - {"the", "a", "an", "in", "of", "to", "and", "or"}
                    pheno_words = set(phenotype.split()) - {"the", "a", "an", "in", "of", "to", "and", "or"}
                    if desc_words and pheno_words:
                        word_overlap = len(desc_words & pheno_words) / max(len(desc_words), 1)
                        confidence += min(0.1, word_overlap * 0.1)

                # Only keep matches with confidence above gene-only baseline
                if confidence > 0.3:
                    scored.append((variant, confidence))

            # Sort by confidence desc, take top N
            scored.sort(key=lambda x: -x[1])
            for variant, confidence in scored[:max_variants_per_entry]:
                match_type = "gene+edit" if confidence >= 0.5 else "gene+partial"
                match = PegRNAClinVarMatch(
                    pegrna_entry_id=entry.id,
                    clinvar_variant_id=variant.id,
                    match_type=match_type,
                    match_confidence=round(confidence, 3),
                )
                session.add(match)
                gene_matches += 1

        if gene_matches > 0:
            session.commit()
            total_matches += gene_matches

        if (i + 1) % 100 == 0:
            console.print(f"  Processed {i+1}/{len(common_genes)} genes, {total_matches:,} matches so far...")

    console.print(f"[green]Created {total_matches:,} ClinVar matches across {len(common_genes)} genes[/green]")
    return total_matches
