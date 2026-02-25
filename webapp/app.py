"""Streamlit web app for the pegRNA database."""
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

import config
from database.models import init_db, Paper, PegRNAEntry
from database.operations import search_entries, get_stats

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

    col5, col6, col7, col8 = st.columns(4)
    with col5:
        pegrna_type_filter = st.selectbox("pegRNA Type", ["All", "pegRNA", "epegRNA"])
    with col6:
        min_eff = st.number_input("Min Efficiency (%)", min_value=0.0, max_value=100.0, value=0.0)
    with col7:
        max_eff = st.number_input("Max Efficiency (%)", min_value=0.0, max_value=100.0, value=100.0)
    with col8:
        func_effect_options = [
            "All",
            "Nonsense (Stop Codon *)",
            "Frameshift",
            "Splice disruption",
            "Knockout",
            "Any LoF",
        ]
        func_effect_filter = st.selectbox("Functional Effect", func_effect_options)

    col9, col10 = st.columns(2)
    with col9:
        validated_only = st.checkbox("Validated entries only")
    with col10:
        detection_options = ["All", "Paper-annotated", "RTT-predicted (computational)"]
        detection_filter = st.selectbox("Stop Codon Detection Method", detection_options)

    # Map UI labels to DB values
    _effect_map = {
        "Nonsense (Stop Codon *)": "Nonsense",
        "Frameshift": "Frameshift",
        "Splice disruption": "Splice disruption",
        "Knockout": "Knockout",
        "Any LoF": "Any LoF",
    }
    func_effect_db = _effect_map.get(func_effect_filter)

    # Query
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
        functional_effect=func_effect_db,
        limit=500,
    )

    # Post-filter by detection method if needed
    if detection_filter == "Paper-annotated":
        results = [
            e for e in results
            if e.functional_effect_detail and not e.functional_effect_detail.startswith("RTT-predicted")
        ]
    elif detection_filter == "RTT-predicted (computational)":
        results = [
            e for e in results
            if e.functional_effect_detail and e.functional_effect_detail.startswith("RTT-predicted")
        ]

    st.markdown(f"**{len(results)} entries found**")

    if results:
        # Build display dataframe
        rows = []
        for e in results:
            # Determine stop codon display
            if e.functional_effect == "Nonsense":
                stop_codon_display = e.functional_effect_detail or "Yes"
                if e.functional_effect_detail and e.functional_effect_detail.startswith("RTT-predicted"):
                    detection_display = "Computational"
                else:
                    detection_display = "Paper-annotated"
            else:
                stop_codon_display = "-"
                detection_display = "-"

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
                "Effect": e.functional_effect or "-",
                "Stop Codon *": stop_codon_display,
                "Detection": detection_display,
                "PE": e.prime_editor or "-",
                "Cell": e.cell_type or "-",
                "Efficiency (%)": f"{float(e.editing_efficiency):.1f}" if e.editing_efficiency is not None else "-",
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
                        if entry.target_region:
                            st.write(f"Target Region: {entry.target_region}"
                                     f"{' - ' + entry.target_region_detail if entry.target_region_detail else ''}")
                        if entry.functional_effect:
                            effect_text = f"**Functional Effect: {entry.functional_effect}**"
                            if entry.functional_effect_detail:
                                effect_text += f" ({entry.functional_effect_detail})"
                            if entry.functional_effect == "Nonsense":
                                if entry.functional_effect_detail and entry.functional_effect_detail.startswith("RTT-predicted"):
                                    effect_text += " — *Computationally detected*"
                                else:
                                    effect_text += " — *Paper-annotated*"
                            st.markdown(effect_text)
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

    # Charts
    if stats["total_entries"] > 0:
        entries = session.query(PegRNAEntry).all()
        df = pd.DataFrame([{
            "edit_type": e.edit_type,
            "pegrna_type": e.pegrna_type,
            "prime_editor": e.prime_editor,
            "cell_type": e.cell_type,
            "target_gene": e.target_gene,
            "efficiency": e.editing_efficiency,
            "year": e.paper.year if e.paper else None,
        } for e in entries])

        st.subheader("Distribution of Edit Types")
        if "edit_type" in df.columns and df["edit_type"].notna().any():
            fig = px.pie(
                df[df["edit_type"].notna()],
                names="edit_type",
                title="Edit Types",
            )
            st.plotly_chart(fig, use_container_width=True)

        st.subheader("Editing Efficiency Distribution")
        eff_data = df[df["efficiency"].notna()]
        if not eff_data.empty:
            fig = px.histogram(
                eff_data,
                x="efficiency",
                nbins=50,
                title="Editing Efficiency Distribution (%)",
                color="pegrna_type",
                barmode="overlay",
            )
            st.plotly_chart(fig, use_container_width=True)

        st.subheader("Top Target Genes")
        if df["target_gene"].notna().any():
            gene_counts = df["target_gene"].value_counts().head(20)
            fig = px.bar(
                x=gene_counts.index,
                y=gene_counts.values,
                title="Top 20 Target Genes",
                labels={"x": "Gene", "y": "Count"},
            )
            st.plotly_chart(fig, use_container_width=True)

        st.subheader("Prime Editor Usage")
        if df["prime_editor"].notna().any():
            pe_counts = df["prime_editor"].value_counts().head(10)
            fig = px.bar(
                x=pe_counts.index,
                y=pe_counts.values,
                title="Prime Editor Usage",
                labels={"x": "Editor", "y": "Count"},
            )
            st.plotly_chart(fig, use_container_width=True)

        # pegRNA vs epegRNA efficiency comparison
        st.subheader("pegRNA vs epegRNA Efficiency")
        type_eff = df[df["efficiency"].notna() & df["pegrna_type"].notna()]
        if not type_eff.empty:
            fig = px.box(
                type_eff,
                x="pegrna_type",
                y="efficiency",
                title="Editing Efficiency: pegRNA vs epegRNA",
                color="pegrna_type",
            )
            st.plotly_chart(fig, use_container_width=True)


# ============================================================
# SUBMIT PAPER PAGE
# ============================================================
elif page == "Submit Paper":
    st.title("Submit Papers for Processing")
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

    if st.button("Process Papers", type="primary"):
        if not input_text.strip():
            st.warning("Please enter at least one identifier.")
        else:
            identifiers = [line.strip() for line in input_text.splitlines() if line.strip()]

            with st.spinner(f"Processing {len(identifiers)} papers..."):
                from discovery.manual_input import parse_identifier, resolve_to_pmid
                from discovery.pubmed_search import fetch_paper_metadata
                from database.operations import get_or_create_paper, update_paper_status

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
                            st.warning(f"Could not resolve: {ident}")
                    else:
                        st.warning(f"Unrecognized: {ident}")

                if pmids:
                    papers = fetch_paper_metadata(pmids)
                    for paper_data in papers:
                        paper, created = get_or_create_paper(session, **paper_data)
                        session.commit()
                        status = "new" if created else "existing"
                        st.write(f"[{status}] {paper.title}")

                    st.success(f"Added {len(papers)} papers. Run the CLI to process them:")
                    st.code("python cli.py extract")
                else:
                    st.error("No valid identifiers found.")

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

This tool was built to help the prime editing community quickly find which pegRNAs
work for their target of interest. Contributions and feedback welcome!
    """)
