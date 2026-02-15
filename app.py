"""Streamlit web app for the pegRNA database - HF Spaces Entry Point."""
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

import re
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

import config
from database.models import init_db, Paper, PegRNAEntry
from database.operations import search_entries, get_stats, sequence_search

# --- Page Config ---
st.set_page_config(
    page_title="pegRNA Database",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- Database Connection ---
@st.cache_resource
def get_session():
    Session = init_db(str(config.DATABASE_PATH))
    return Session()


session = get_session()

# --- Sidebar Navigation ---
st.sidebar.title("pegRNA Database")
st.sidebar.markdown("*Experimentally validated (e)pegRNAs from the literature*")

page = st.sidebar.radio(
    "Navigate",
    ["Search & Browse", "Statistics", "Submit Paper", "Export", "About"],
)


# ============================================================
# SEARCH & BROWSE PAGE
# ============================================================
if page == "Search & Browse":
    st.title("Search (e)pegRNA Entries")

    # Filters
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        gene_filter = st.text_input("Target Gene", placeholder="e.g., HEK3, FANCF")
    with col2:
        edit_type_filter = st.selectbox(
            "Edit Type",
            ["All", "substitution", "insertion", "deletion"],
        )
    with col3:
        pe_filter = st.text_input("Prime Editor", placeholder="e.g., PE2, PE3, PE7")
    with col4:
        cell_filter = st.text_input("Cell Type", placeholder="e.g., HEK293T, K562")

    col5, col6, col7 = st.columns(3)
    with col5:
        pegrna_type_filter = st.selectbox("pegRNA Type", ["All", "pegRNA", "epegRNA"])
    with col6:
        min_eff = st.number_input("Min Efficiency (%)", min_value=0.0, max_value=100.0, value=0.0)
    with col7:
        max_eff = st.number_input("Max Efficiency (%)", min_value=0.0, max_value=100.0, value=100.0)

    validated_only = st.checkbox("Validated entries only")

    # --- Sequence Search ---
    with st.expander("Sequence Search (find similar pegRNAs by DNA sequence)", expanded=False):
        seq_col1, seq_col2, seq_col3 = st.columns([3, 1, 1])
        with seq_col1:
            query_seq = st.text_input(
                "Query Sequence",
                placeholder="e.g., GGCCCAGACTGAGCACGTGA",
                help="Enter a DNA sequence (ACGT) to find similar pegRNAs",
            )
        with seq_col2:
            seq_field = st.selectbox("Search Field", ["spacer", "pbs", "rtt", "full"])
        with seq_col3:
            max_dist = st.selectbox("Max Mismatches", [0, 1, 2, 3, 4, 5], index=2)

        if query_seq:
            clean_seq = query_seq.strip().upper().replace("U", "T").replace(" ", "")
            if not re.match(r"^[ACGT]+$", clean_seq):
                st.error("Invalid sequence. Use only A, C, G, T characters.")
            elif len(clean_seq) < 10:
                st.warning("Sequence too short. Enter at least 10 nucleotides.")
            else:
                with st.spinner(f"Searching {seq_field} sequences..."):
                    seq_results = sequence_search(
                        session, clean_seq,
                        max_distance=max_dist,
                        search_field=seq_field,
                        limit=200,
                    )
                if seq_results:
                    st.success(f"Found {len(seq_results)} matches (within {max_dist} mismatches)")
                    seq_rows = []
                    for r in seq_results:
                        e = r["entry"]
                        attr = f"{seq_field}_sequence" if seq_field != "full" else "full_sequence"
                        seq_rows.append({
                            "Dist": r["distance"],
                            "ID": e.id,
                            "Name": e.entry_name or "-",
                            "Gene": e.target_gene or "-",
                            f"{seq_field.upper()}": getattr(e, attr) or "-",
                            "Edit": e.edit_type or "-",
                            "Eff (%)": f"{e.editing_efficiency:.1f}" if e.editing_efficiency else "-",
                            "PE": e.prime_editor or "-",
                            "Cell": e.cell_type or "-",
                            "PMID": e.paper.pmid if e.paper else "-",
                        })
                    st.dataframe(pd.DataFrame(seq_rows), use_container_width=True, height=400)
                else:
                    st.info("No matches found. Try increasing the max mismatches.")

    col_page1, col_page2 = st.columns(2)
    with col_page1:
        page_size = st.selectbox("Results per page", [100, 500, 1000, 5000], index=1)
    with col_page2:
        page_num = st.number_input("Page", min_value=1, value=1, step=1)

    # Count total matching entries
    from database.models import PegRNAEntry as _PE
    count_query = session.query(_PE)
    if gene_filter:
        count_query = count_query.filter(_PE.target_gene.ilike(f"%{gene_filter}%"))
    if edit_type_filter != "All":
        count_query = count_query.filter(_PE.edit_type.ilike(f"%{edit_type_filter}%"))
    if pe_filter:
        count_query = count_query.filter(_PE.prime_editor.ilike(f"%{pe_filter}%"))
    if cell_filter:
        count_query = count_query.filter(_PE.cell_type.ilike(f"%{cell_filter}%"))
    if pegrna_type_filter != "All":
        count_query = count_query.filter(_PE.pegrna_type.ilike(f"%{pegrna_type_filter}%"))
    if min_eff > 0:
        count_query = count_query.filter(_PE.editing_efficiency >= min_eff)
    if max_eff < 100:
        count_query = count_query.filter(_PE.editing_efficiency <= max_eff)
    if validated_only:
        count_query = count_query.filter(_PE.validated.is_(True))
    total_count = count_query.count()
    total_pages = max(1, (total_count + page_size - 1) // page_size)

    # Query
    offset = (page_num - 1) * page_size
    results = search_entries(
        session,
        target_gene=gene_filter or None,
        edit_type=edit_type_filter if edit_type_filter != "All" else None,
        prime_editor=pe_filter or None,
        cell_type=cell_filter or None,
        pegrna_type=pegrna_type_filter if pegrna_type_filter != "All" else None,
        min_efficiency=min_eff if min_eff > 0 else None,
        max_efficiency=max_eff if max_eff < 100 else None,
        validated_only=validated_only,
        limit=page_size,
        offset=offset,
    )

    st.markdown(f"**{total_count:,} entries found** (showing {offset+1}-{min(offset+page_size, total_count)} of {total_count:,}, page {page_num}/{total_pages})")

    if results:
        # Build display dataframe
        rows = []
        for e in results:
            rows.append({
                "ID": e.id,
                "Name": e.entry_name or "-",
                "Type": e.pegrna_type or "-",
                "Gene": e.target_gene or "-",
                "Spacer": e.spacer_sequence or "-",
                "PBS": f"{e.pbs_length}nt" if e.pbs_length else "-",
                "RTT": f"{e.rtt_length}nt" if e.rtt_length else "-",
                "3' Motif": e.three_prime_extension or "-",
                "Edit": e.edit_type or "-",
                "PE": e.prime_editor or "-",
                "Cell": e.cell_type or "-",
                "Efficiency (%)": f"{e.editing_efficiency:.1f}" if e.editing_efficiency is not None else "-",
                "Confidence": f"{e.confidence_score:.2f}" if e.confidence_score else "-",
                "Paper PMID": e.paper.pmid if e.paper else "-",
            })

        df = pd.DataFrame(rows)
        st.dataframe(df, use_container_width=True, height=600)

        # Detail expander for selected entry
        selected_id = st.selectbox("View entry details", [r["ID"] for r in rows])
        if selected_id:
            entry = session.query(PegRNAEntry).get(selected_id)
            if entry:
                with st.expander("Entry Details", expanded=True):
                    dcol1, dcol2 = st.columns(2)
                    with dcol1:
                        st.markdown("**Sequences**")
                        st.code(f"Spacer:  {entry.spacer_sequence or 'N/A'}")
                        st.code(f"PBS:     {entry.pbs_sequence or 'N/A'} ({entry.pbs_length or '?'}nt)")
                        st.code(f"RTT:     {entry.rtt_sequence or 'N/A'} ({entry.rtt_length or '?'}nt)")
                        if entry.three_prime_extension:
                            st.code(f"3' ext:  {entry.three_prime_extension}")
                        if entry.full_sequence:
                            st.code(f"Full:    {entry.full_sequence}")
                        if entry.nicking_sgrna_seq:
                            st.code(f"Nick sgRNA: {entry.nicking_sgrna_seq}")

                    with dcol2:
                        st.markdown("**Experimental Details**")
                        st.write(f"Target: {entry.target_gene} ({entry.target_organism or 'N/A'})")
                        st.write(f"Edit: {entry.edit_type} - {entry.edit_description or 'N/A'}")
                        st.write(f"Prime Editor: {entry.prime_editor or 'N/A'}")
                        st.write(f"Cell Type: {entry.cell_type or 'N/A'}")
                        st.write(f"Delivery: {entry.delivery_method or 'N/A'}")
                        st.markdown("**Results**")
                        st.write(f"Editing Efficiency: {entry.editing_efficiency or 'N/A'}%")
                        st.write(f"Product Purity: {entry.product_purity or 'N/A'}%")
                        st.write(f"Indel Frequency: {entry.indel_frequency or 'N/A'}%")

                    if entry.paper:
                        st.markdown("**Source**")
                        st.write(f"{entry.paper.title}")
                        st.write(f"{entry.paper.authors}")
                        st.write(f"{entry.paper.journal} ({entry.paper.year})")
                        if entry.paper.pmid:
                            st.markdown(
                                f"[PubMed](https://pubmed.ncbi.nlm.nih.gov/{entry.paper.pmid}/)"
                            )
    else:
        st.info("No entries match your filters. Try broadening your search.")


# ============================================================
# STATISTICS PAGE
# ============================================================
elif page == "Statistics":
    st.title("Database Statistics")

    stats = get_stats(session)

    # KPI cards
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total Papers", stats["total_papers"])
    col2.metric("Processed Papers", stats["completed_papers"])
    col3.metric("Total Entries", stats["total_entries"])
    col4.metric("Validated Entries", stats["validated_entries"])

    col5, col6, col7 = st.columns(3)
    col5.metric("Unique Genes", stats["unique_genes"])
    col6.metric("Unique Cell Types", stats["unique_cell_types"])
    col7.metric("Unique Organisms", stats["unique_organisms"])

    # Charts - using SQL aggregations for performance (264k+ entries)
    if stats["total_entries"] > 0:
        from sqlalchemy import func, distinct, case

        # Edit type distribution
        st.subheader("Distribution of Edit Types")
        edit_type_counts = (
            session.query(PegRNAEntry.edit_type, func.count(PegRNAEntry.id))
            .filter(PegRNAEntry.edit_type.isnot(None))
            .group_by(PegRNAEntry.edit_type)
            .all()
        )
        if edit_type_counts:
            et_df = pd.DataFrame(edit_type_counts, columns=["edit_type", "count"])
            fig = px.pie(et_df, names="edit_type", values="count", title="Edit Types")
            st.plotly_chart(fig, use_container_width=True)

        # Efficiency distribution - sample for histogram
        st.subheader("Editing Efficiency Distribution")
        eff_rows = (
            session.query(PegRNAEntry.editing_efficiency, PegRNAEntry.pegrna_type)
            .filter(PegRNAEntry.editing_efficiency.isnot(None))
            .limit(50000)
            .all()
        )
        if eff_rows:
            eff_df = pd.DataFrame(eff_rows, columns=["efficiency", "pegrna_type"])
            fig = px.histogram(
                eff_df, x="efficiency", nbins=50,
                title=f"Editing Efficiency Distribution (%) - {len(eff_rows):,} entries",
                color="pegrna_type", barmode="overlay",
            )
            st.plotly_chart(fig, use_container_width=True)

        # Top genes
        st.subheader("Top Target Genes")
        gene_counts = (
            session.query(PegRNAEntry.target_gene, func.count(PegRNAEntry.id).label("count"))
            .filter(PegRNAEntry.target_gene.isnot(None))
            .group_by(PegRNAEntry.target_gene)
            .order_by(func.count(PegRNAEntry.id).desc())
            .limit(20)
            .all()
        )
        if gene_counts:
            gene_df = pd.DataFrame(gene_counts, columns=["Gene", "Count"])
            fig = px.bar(gene_df, x="Gene", y="Count", title="Top 20 Target Genes")
            st.plotly_chart(fig, use_container_width=True)

        # Prime editor usage
        st.subheader("Prime Editor Usage")
        pe_counts = (
            session.query(PegRNAEntry.prime_editor, func.count(PegRNAEntry.id).label("count"))
            .filter(PegRNAEntry.prime_editor.isnot(None))
            .group_by(PegRNAEntry.prime_editor)
            .order_by(func.count(PegRNAEntry.id).desc())
            .limit(10)
            .all()
        )
        if pe_counts:
            pe_df = pd.DataFrame(pe_counts, columns=["Editor", "Count"])
            fig = px.bar(pe_df, x="Editor", y="Count", title="Prime Editor Usage")
            st.plotly_chart(fig, use_container_width=True)

        # pegRNA vs epegRNA efficiency - use aggregated stats
        st.subheader("pegRNA vs epegRNA Efficiency")
        type_eff_rows = (
            session.query(PegRNAEntry.editing_efficiency, PegRNAEntry.pegrna_type)
            .filter(
                PegRNAEntry.editing_efficiency.isnot(None),
                PegRNAEntry.pegrna_type.isnot(None),
            )
            .limit(50000)
            .all()
        )
        if type_eff_rows:
            type_eff_df = pd.DataFrame(type_eff_rows, columns=["efficiency", "pegrna_type"])
            fig = px.box(
                type_eff_df, x="pegrna_type", y="efficiency",
                title="Editing Efficiency: pegRNA vs epegRNA",
                color="pegrna_type",
            )
            st.plotly_chart(fig, use_container_width=True)

        # Papers by year
        st.subheader("Papers by Year")
        year_counts = (
            session.query(Paper.year, func.count(Paper.id).label("count"))
            .filter(Paper.year.isnot(None), Paper.extraction_status == "completed")
            .group_by(Paper.year)
            .order_by(Paper.year)
            .all()
        )
        if year_counts:
            year_df = pd.DataFrame(year_counts, columns=["Year", "Papers"])
            fig = px.bar(year_df, x="Year", y="Papers", title="Completed Papers by Year")
            st.plotly_chart(fig, use_container_width=True)


# ============================================================
# SUBMIT PAPER PAGE
# ============================================================
elif page == "Submit Paper":
    st.title("Submit Papers for Processing")

    has_ncbi = bool(config.NCBI_EMAIL)
    has_llm = bool(config.ANTHROPIC_API_KEY or config.OPENROUTER_API_KEY or config.CHATGPT_API_KEY)

    if not has_ncbi:
        st.warning(
            "NCBI_EMAIL is not configured. Paper retrieval from PubMed/PMC requires "
            "an email address. Set it as a secret in your HF Space settings or in `.env`."
        )

    st.markdown(
        "Enter PMIDs, DOIs, or PubMed/Nature/PMC URLs to extract pegRNA data. "
        "One identifier per line."
    )

    input_text = st.text_area(
        "Paper Identifiers",
        height=150,
        placeholder="34608327\nhttps://www.nature.com/articles/s41586-024-07259-6\n10.1038/s41587-022-01613-7",
    )

    method = st.selectbox("Extraction Method", ["auto", "rule", "llm"])

    if not has_llm and method in ("auto", "llm"):
        st.info("No LLM API key configured. LLM fallback extraction will be skipped. "
                "Rule-based extraction from supplementary tables will still work.")

    if st.button("Process Papers", type="primary"):
        if not input_text.strip():
            st.warning("Please enter at least one identifier.")
        else:
            identifiers = [line.strip() for line in input_text.splitlines() if line.strip()]

            progress = st.empty()
            log_area = st.container()

            try:
                from discovery.manual_input import parse_identifier, resolve_to_pmid
                from discovery.pubmed_search import fetch_paper_metadata
                from database.operations import get_or_create_paper, update_paper_status, bulk_add_entries
                from database.models import Paper

                # Step 1: Resolve identifiers to PMIDs
                progress.info("Resolving identifiers...")
                pmids = []
                for ident in identifiers:
                    parsed = parse_identifier(ident)
                    if parsed["type"] == "pmid":
                        pmids.append(parsed["value"])
                    elif parsed["type"] in ("doi", "pmcid"):
                        pmid = resolve_to_pmid(parsed)
                        if pmid:
                            pmids.append(pmid)
                        else:
                            log_area.warning(f"Could not resolve: {ident}")
                    else:
                        log_area.warning(f"Unrecognized: {ident}")

                if not pmids:
                    st.error("No valid identifiers found.")
                    st.stop()

                # Step 2: Fetch metadata and add to database
                progress.info("Fetching paper metadata from PubMed...")
                papers_data = fetch_paper_metadata(pmids)
                papers_to_process = []

                for paper_data in papers_data:
                    paper, created = get_or_create_paper(session, **paper_data)
                    session.commit()
                    status_label = "New" if created else "Existing"
                    log_area.write(f"**[{status_label}]** {paper.title}")
                    if created or paper.extraction_status == "pending":
                        papers_to_process.append(paper)

                if not papers_to_process:
                    st.success("All papers already processed! Browse them in Search & Browse.")
                    st.stop()

                # Step 3: Process each paper through the extraction pipeline
                from retrieval.pmc_fetcher import fetch_full_text, extract_text_from_bioc, fetch_by_pmid
                from retrieval.supplementary import (
                    list_supplementary_files, download_supplementary_file,
                    parse_supplementary_tables, find_pegrna_tables,
                )
                from extraction.rule_based import extract_from_multiple_tables, extract_from_text
                from extraction.validator import validate_batch

                total_entries = 0

                for i, paper in enumerate(papers_to_process):
                    progress.info(f"Processing paper {i+1}/{len(papers_to_process)}: {paper.pmid}...")
                    all_entries = []
                    paper_text = ""

                    # Fetch full text
                    bioc_data = None
                    if paper.pmcid:
                        bioc_data = fetch_full_text(paper.pmcid)
                    elif paper.pmid:
                        bioc_data = fetch_by_pmid(paper.pmid)
                    if bioc_data:
                        paper_text = extract_text_from_bioc(bioc_data)

                    # Try supplementary materials (rule-based)
                    if paper.pmcid and method in ("auto", "rule"):
                        log_area.write(f"Fetching supplementary files for {paper.pmcid}...")
                        supp_files = list_supplementary_files(paper.pmcid)
                        tabular_files = [f for f in supp_files if f["type"] in ("excel", "csv", "tsv")]

                        for supp in tabular_files:
                            filepath = download_supplementary_file(paper.pmcid, supp["filename"])
                            if filepath:
                                dfs = parse_supplementary_tables(filepath)
                                pegrna_tables = find_pegrna_tables(dfs)
                                if pegrna_tables:
                                    entries = extract_from_multiple_tables(pegrna_tables)
                                    all_entries.extend(entries)
                                    log_area.write(f"Found {len(entries)} entries in {supp['filename']}")

                    # Text-based extraction
                    if paper_text and method in ("auto", "rule"):
                        text_entries = extract_from_text(paper_text)
                        all_entries.extend(text_entries)

                    # LLM fallback
                    if not all_entries and paper_text and has_llm and method in ("auto", "llm"):
                        log_area.write("Rule-based found nothing, trying LLM extraction...")
                        try:
                            from extraction.llm_extractor import extract_with_llm
                            llm_entries = extract_with_llm(
                                paper_text,
                                paper_title=paper.title or "",
                                paper_abstract=paper.abstract or "",
                            )
                            all_entries.extend(llm_entries)
                        except Exception as llm_err:
                            log_area.warning(f"LLM extraction failed: {llm_err}")

                    # Validate and store
                    if all_entries:
                        valid_entries, invalid = validate_batch(all_entries)
                        if valid_entries:
                            entry_dicts = [e.to_db_dict() for e in valid_entries]
                            count = bulk_add_entries(session, paper.id, entry_dicts)
                            update_paper_status(session, paper.id, "completed", method=method)
                            session.commit()
                            total_entries += count
                            log_area.success(f"Stored **{count}** entries for {paper.pmid}")
                        else:
                            update_paper_status(session, paper.id, "review",
                                notes=f"All {len(all_entries)} entries failed validation")
                            session.commit()
                            log_area.warning(f"All entries failed validation for {paper.pmid}")
                    else:
                        update_paper_status(session, paper.id, "failed",
                            notes="No pegRNA data could be extracted")
                        session.commit()
                        log_area.error(f"No pegRNA data found in {paper.pmid}")

                progress.empty()
                if total_entries > 0:
                    st.success(f"Done! Extracted **{total_entries}** pegRNA entries from "
                              f"{len(papers_to_process)} papers. Browse them in Search & Browse.")
                else:
                    st.warning("Processing complete but no pegRNA entries were extracted. "
                              "The papers may not contain structured pegRNA data in their "
                              "supplementary materials.")

            except Exception as e:
                progress.empty()
                st.error(f"Error: {str(e)}")
                import traceback
                st.code(traceback.format_exc())

    # Show papers in database
    st.subheader("Papers in Database")
    papers = session.query(Paper).order_by(Paper.year.desc()).all()
    if papers:
        paper_data = [{
            "PMID": p.pmid or "-",
            "Title": (p.title or "")[:80],
            "Year": p.year,
            "Journal": (p.journal or "")[:30],
            "Status": p.extraction_status,
            "Entries": len(p.entries),
        } for p in papers]
        st.dataframe(pd.DataFrame(paper_data), use_container_width=True)
    else:
        st.info("No papers in database yet.")


# ============================================================
# EXPORT PAGE
# ============================================================
elif page == "Export":
    st.title("Export Database")

    export_format = st.selectbox("Format", ["CSV", "JSON"])
    validated_only = st.checkbox("Validated entries only", value=False)

    entries = session.query(PegRNAEntry).all()
    if validated_only:
        entries = [e for e in entries if e.validated]

    st.write(f"**{len(entries)} entries** to export")

    if entries and st.button("Generate Export"):
        rows = []
        for e in entries:
            rows.append({
                "entry_name": e.entry_name,
                "pegrna_type": e.pegrna_type,
                "spacer_sequence": e.spacer_sequence,
                "pbs_sequence": e.pbs_sequence,
                "pbs_length": e.pbs_length,
                "rtt_sequence": e.rtt_sequence,
                "rtt_length": e.rtt_length,
                "three_prime_extension": e.three_prime_extension,
                "full_sequence": e.full_sequence,
                "nicking_sgrna_seq": e.nicking_sgrna_seq,
                "target_gene": e.target_gene,
                "target_locus": e.target_locus,
                "target_organism": e.target_organism,
                "edit_type": e.edit_type,
                "edit_description": e.edit_description,
                "intended_mutation": e.intended_mutation,
                "prime_editor": e.prime_editor,
                "cell_type": e.cell_type,
                "delivery_method": e.delivery_method,
                "editing_efficiency": e.editing_efficiency,
                "product_purity": e.product_purity,
                "indel_frequency": e.indel_frequency,
                "confidence_score": e.confidence_score,
                "validated": e.validated,
                "paper_pmid": e.paper.pmid if e.paper else None,
                "paper_doi": e.paper.doi if e.paper else None,
                "paper_year": e.paper.year if e.paper else None,
            })

        df = pd.DataFrame(rows)

        if export_format == "CSV":
            csv_data = df.to_csv(index=False)
            st.download_button(
                "Download CSV",
                csv_data,
                file_name="pegrna_database.csv",
                mime="text/csv",
            )
        else:
            json_data = df.to_json(orient="records", indent=2)
            st.download_button(
                "Download JSON",
                json_data,
                file_name="pegrna_database.json",
                mime="application/json",
            )


# ============================================================
# ABOUT PAGE
# ============================================================
elif page == "About":
    st.title("About the pegRNA Database")
    st.markdown("""
## What is this?

This is an automated database of experimentally validated prime editing guide RNAs
(pegRNAs and epegRNAs) extracted from the scientific literature.

**Prime editing** is a CRISPR-based genome editing technology that uses a
prime editing guide RNA (pegRNA) to direct precise edits without requiring
double-strand breaks or donor DNA templates.

## How does it work?

1. **Paper Discovery**: We search PubMed for papers containing prime editing experiments
2. **Data Retrieval**: Full text and supplementary materials are fetched via the PMC BioC API
3. **Extraction**: A hybrid pipeline uses rule-based parsing (for structured tables) and
   LLM-powered extraction (for unstructured text) to identify pegRNA data
4. **Validation**: Extracted data is validated (sequence checks, range checks)
5. **Database**: All validated entries are stored and made searchable

## Data Fields

Each entry contains:
- **Sequences**: spacer, PBS, RTT, 3' extension, full sequence, nicking sgRNA
- **Target**: gene, genomic locus, organism
- **Edit**: type (substitution/insertion/deletion), description, intended mutation
- **Conditions**: prime editor version, cell type, delivery method
- **Results**: editing efficiency, product purity, indel frequency
- **Provenance**: source paper, confidence score, validation status

## Confidence Scoring

- **0.8+**: Rule-based extraction from structured supplementary tables (highest confidence)
- **0.5-0.7**: LLM-based extraction from paper text (moderate confidence)
- **Validated**: Human-reviewed entries flagged as validated

## Source Code

This database is available on [GitHub](https://github.com/JoeyRiepsaame/pegrna-db) and
hosted on [Hugging Face Spaces](https://huggingface.co/spaces/JoeyRiepsaame/pegRNA-DB).

Built to help the prime editing community quickly find which pegRNAs work for their target of interest.
Contributions and feedback welcome!

## Citation

If you use this database in your research, please cite:
```
pegRNA Database (2026)
https://huggingface.co/spaces/JoeyRiepsaame/pegRNA-DB
```
    """)
