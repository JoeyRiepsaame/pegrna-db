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
import numpy as np

import config
from database.models import init_db, Paper, PegRNAEntry
from database.operations import search_entries, get_stats, sequence_search

def strip_html(text: str) -> str:
    """Strip HTML tags from text (e.g., <i>in vivo</i> -> in vivo)."""
    return re.sub(r"<[^>]+>", "", text) if text else text

# SpCas9 scaffold (tracrRNA) used for full pegRNA sequence reconstruction
SCAFFOLD = "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"

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


@st.cache_data(ttl=300)
def get_distinct_pe_versions():
    """Get distinct prime editor values from DB for dropdown."""
    from sqlalchemy import distinct
    s = get_session()
    values = s.query(distinct(PegRNAEntry.prime_editor)).filter(
        PegRNAEntry.prime_editor.isnot(None)
    ).all()
    return sorted([v[0] for v in values if v[0]])


def has_clinvar_data(session):
    """Check if ClinVar data has been loaded."""
    try:
        from database.models import ClinVarVariant
        return session.query(ClinVarVariant).limit(1).first() is not None
    except Exception:
        return False


def classify_editing_technology(editor_value):
    """Classify an editor value as Prime Editing or Base Editing."""
    if not editor_value:
        return "Unknown"
    val = editor_value.lower()
    if any(be in val for be in ["abe", "cbe"]):
        return "Base Editing"
    if any(pe in val for pe in ["pe1", "pe2", "pe3", "pe4", "pe5", "pe6", "pe7",
                                 "pemax", "peco", "twipe", "twinpe", "mpe",
                                 "sape", "tj-pe", "epe", "pegr", "prime"]):
        return "Prime Editing"
    return "Unknown"


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
    st.title("Search Guide RNA Entries")

    # Clear filters button
    if st.button("Clear Filters"):
        for key in list(st.session_state.keys()):
            if key.startswith("filter_"):
                del st.session_state[key]
        st.rerun()

    # Filters - Row 1
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        gene_filter = st.text_input("Target Gene", placeholder="e.g., HEK3, FANCF", key="filter_gene")
    with col2:
        edit_type_filter = st.selectbox(
            "Edit Type",
            ["All", "substitution", "insertion", "deletion"],
            key="filter_edit_type",
        )
    with col3:
        pe_options = ["All"] + get_distinct_pe_versions()
        pe_filter = st.selectbox("Prime Editor", pe_options, key="filter_pe")
    with col4:
        cell_filter = st.text_input("Cell Type", placeholder="e.g., HEK293T, K562", key="filter_cell")

    # Filters - Row 2
    col5, col6, col7, col8, col9 = st.columns(5)
    with col5:
        pegrna_type_filter = st.selectbox("pegRNA Type", ["All", "pegRNA", "epegRNA"], key="filter_type")
    with col6:
        tech_filter = st.selectbox("Editing Technology", ["All", "Prime Editing", "Base Editing"], key="filter_tech")
    with col7:
        min_eff = st.number_input("Min Efficiency (%)", min_value=0.0, max_value=100.0, value=0.0, key="filter_min_eff")
    with col8:
        max_eff = st.number_input("Max Efficiency (%)", min_value=0.0, max_value=100.0, value=100.0, key="filter_max_eff")
    with col9:
        organism_filter = st.text_input("Organism", placeholder="e.g., Human, Mouse", key="filter_organism")

    # Filters - Row 3: PMID + Author + Paper Title
    col_pmid, col_author, col_title = st.columns(3)
    with col_pmid:
        pmid_filter = st.text_input("PMID", placeholder="e.g., 36553615", key="filter_pmid")
    with col_author:
        author_filter = st.text_input("Author", placeholder="e.g., Anzalone, Chen", key="filter_author")
    with col_title:
        title_filter = st.text_input("Paper Title", placeholder="e.g., prime editing, high-throughput", key="filter_title")

    # Row 4: validated + ClinVar filter + sorting
    col_v, col_cv, col_sort, col_order = st.columns(4)
    with col_v:
        validated_only = st.checkbox("Validated entries only", key="filter_validated")
    with col_cv:
        clinvar_filter_options = ["All"]
        _has_clinvar = has_clinvar_data(session)
        if _has_clinvar:
            clinvar_filter_options += ["Has ClinVar match", "Pathogenic match", "No ClinVar match"]
        clinvar_filter = st.selectbox("ClinVar Match", clinvar_filter_options, key="filter_clinvar")
    with col_sort:
        sort_by = st.selectbox("Sort by", [None, "Efficiency", "Gene", "Edit Type", "PE Version"], key="filter_sort")
    with col_order:
        sort_desc = st.checkbox("Descending", value=True, key="filter_sort_desc")

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

    # Build count query with all filters
    from database.models import PegRNAEntry as _PE
    count_query = session.query(_PE)
    if pmid_filter or author_filter or title_filter:
        from database.models import Paper as _P
        count_query = count_query.join(_P, _PE.paper_id == _P.id)
    if pmid_filter:
        count_query = count_query.filter(_P.pmid == pmid_filter.strip())
    if author_filter:
        count_query = count_query.filter(_P.authors.ilike(f"%{author_filter}%"))
    if title_filter:
        count_query = count_query.filter(_P.title.ilike(f"%{title_filter}%"))
    if gene_filter:
        count_query = count_query.filter(_PE.target_gene.ilike(f"%{gene_filter}%"))
    if edit_type_filter != "All":
        count_query = count_query.filter(_PE.edit_type.ilike(f"%{edit_type_filter}%"))
    if pe_filter != "All":
        count_query = count_query.filter(_PE.prime_editor.ilike(f"%{pe_filter}%"))
    if cell_filter:
        count_query = count_query.filter(_PE.cell_type.ilike(f"%{cell_filter}%"))
    if pegrna_type_filter != "All":
        count_query = count_query.filter(_PE.pegrna_type.ilike(f"%{pegrna_type_filter}%"))
    if min_eff > 0:
        count_query = count_query.filter(_PE.editing_efficiency >= min_eff)
    if max_eff < 100:
        count_query = count_query.filter(_PE.editing_efficiency <= max_eff)
    if organism_filter:
        count_query = count_query.filter(_PE.target_organism.ilike(f"%{organism_filter}%"))
    if validated_only:
        count_query = count_query.filter(_PE.validated.is_(True))
    if tech_filter == "Prime Editing":
        from sqlalchemy import or_
        count_query = count_query.filter(or_(
            *[_PE.prime_editor.ilike(f"%{kw}%") for kw in
              ["PE1", "PE2", "PE3", "PE4", "PE5", "PE6", "PE7", "PEmax", "PECO", "twinPE", "mPE", "SaPE", "TJ-PE", "ePE", "pegRNA", "prime"]]
        ))
    elif tech_filter == "Base Editing":
        from sqlalchemy import or_
        count_query = count_query.filter(or_(
            _PE.prime_editor.ilike("%ABE%"),
            _PE.prime_editor.ilike("%CBE%"),
        ))

    # ClinVar filter on count query
    if _has_clinvar and clinvar_filter != "All":
        from database.models import PegRNAClinVarMatch, ClinVarVariant
        from sqlalchemy import exists
        match_subq = session.query(PegRNAClinVarMatch.pegrna_entry_id).subquery()
        if clinvar_filter == "Has ClinVar match":
            count_query = count_query.filter(_PE.id.in_(session.query(PegRNAClinVarMatch.pegrna_entry_id)))
        elif clinvar_filter == "Pathogenic match":
            pathogenic_subq = (
                session.query(PegRNAClinVarMatch.pegrna_entry_id)
                .join(ClinVarVariant, PegRNAClinVarMatch.clinvar_variant_id == ClinVarVariant.id)
                .filter(ClinVarVariant.clinical_significance.ilike("%pathogenic%"))
                .subquery()
            )
            count_query = count_query.filter(_PE.id.in_(pathogenic_subq))
        elif clinvar_filter == "No ClinVar match":
            count_query = count_query.filter(~_PE.id.in_(session.query(PegRNAClinVarMatch.pegrna_entry_id)))

    total_count = count_query.count()
    total_pages = max(1, (total_count + page_size - 1) // page_size)

    # Query results
    offset = (page_num - 1) * page_size
    results = search_entries(
        session,
        target_gene=gene_filter or None,
        edit_type=edit_type_filter if edit_type_filter != "All" else None,
        prime_editor=pe_filter if pe_filter != "All" else None,
        cell_type=cell_filter or None,
        pegrna_type=pegrna_type_filter if pegrna_type_filter != "All" else None,
        min_efficiency=min_eff if min_eff > 0 else None,
        max_efficiency=max_eff if max_eff < 100 else None,
        target_organism=organism_filter or None,
        editing_technology=tech_filter if tech_filter != "All" else None,
        validated_only=validated_only,
        pmid=pmid_filter or None,
        author=author_filter or None,
        paper_title=title_filter or None,
        sort_by=sort_by,
        sort_desc=sort_desc,
        limit=page_size,
        offset=offset,
    )

    st.markdown(f"**{total_count:,} entries found** (showing {offset+1}-{min(offset+page_size, total_count)} of {total_count:,}, page {page_num}/{total_pages})")

    if results:
        # Pre-fetch ClinVar match counts for displayed entries
        clinvar_counts = {}
        if _has_clinvar:
            from database.models import PegRNAClinVarMatch
            from sqlalchemy import func
            entry_ids = [e.id for e in results]
            counts = (
                session.query(PegRNAClinVarMatch.pegrna_entry_id, func.count(PegRNAClinVarMatch.id))
                .filter(PegRNAClinVarMatch.pegrna_entry_id.in_(entry_ids))
                .group_by(PegRNAClinVarMatch.pegrna_entry_id)
                .all()
            )
            clinvar_counts = dict(counts)

        # Build display dataframe
        rows = []
        for e in results:
            ext = getattr(e, 'extension_sequence', None)
            if not ext and e.rtt_sequence and e.pbs_sequence:
                ext = e.rtt_sequence + e.pbs_sequence
            tech = classify_editing_technology(e.prime_editor)
            row = {
                "ID": e.id,
                "Name": e.entry_name or "-",
                "Tech": tech,
                "Type": e.pegrna_type or "-",
                "Gene": e.target_gene or "-",
                "Spacer": e.spacer_sequence or "-",
                "3' Extension (RTT + PBS)": ext or "-",
                "PBS": f"{e.pbs_length}nt" if e.pbs_length else "-",
                "RTT": f"{e.rtt_length}nt" if e.rtt_length else "-",
                "Edit": e.edit_type or "-",
                "Editor": e.prime_editor or "-",
                "Organism": e.target_organism or "-",
                "Efficiency (%)": f"{e.editing_efficiency:.1f}" if e.editing_efficiency is not None else "-",
                "Paper PMID": e.paper.pmid if e.paper else "-",
                "Paper Title": (strip_html(e.paper.title)[:60] + "...") if e.paper and e.paper.title and len(strip_html(e.paper.title)) > 60 else (strip_html(e.paper.title) if e.paper and e.paper.title else "-"),
            }
            if _has_clinvar:
                row["ClinVar"] = clinvar_counts.get(e.id, 0)
            rows.append(row)

        df = pd.DataFrame(rows)
        st.dataframe(df, use_container_width=True, height=600)

        # --- Export filtered results ---
        st.markdown("**Export these results**")
        exp_col1, exp_col2, exp_col3 = st.columns(3)

        with exp_col1:
            csv_export = df.to_csv(index=False)
            st.download_button(
                f"Download CSV ({len(rows)} rows)",
                csv_export,
                file_name="pegrna_search_results.csv",
                mime="text/csv",
            )

        with exp_col2:
            # FASTA export
            fasta_lines = []
            for e in results:
                ext = getattr(e, 'extension_sequence', None)
                if not ext and e.rtt_sequence and e.pbs_sequence:
                    ext = e.rtt_sequence + e.pbs_sequence
                full = e.full_sequence
                if not full and e.spacer_sequence and ext:
                    full = e.spacer_sequence + SCAFFOLD + ext
                if full:
                    header = f">{e.id}|{e.target_gene or 'unknown'}|{e.edit_type or 'unknown'}|eff={e.editing_efficiency or 'NA'}%"
                    fasta_lines.append(header)
                    fasta_lines.append(full)
            if fasta_lines:
                st.download_button(
                    f"Download FASTA ({len(fasta_lines)//2} seqs)",
                    "\n".join(fasta_lines),
                    file_name="pegrna_search_results.fasta",
                    mime="text/plain",
                )
            else:
                st.caption("No full sequences available for FASTA export")

        with exp_col3:
            # SnapGene GenBank format
            gb_records = []
            for e in results:
                ext = getattr(e, 'extension_sequence', None)
                if not ext and e.rtt_sequence and e.pbs_sequence:
                    ext = e.rtt_sequence + e.pbs_sequence
                full = e.full_sequence
                if not full and e.spacer_sequence and ext:
                    full = e.spacer_sequence + SCAFFOLD + ext
                if not full:
                    continue

                locus_name = f"pegRNA_{e.id}"
                gene = e.target_gene or "unknown"
                seq_len = len(full)

                # Find component positions
                spacer_end = len(e.spacer_sequence) if e.spacer_sequence else 20
                scaffold_end = spacer_end + len(SCAFFOLD)
                rtt_end = scaffold_end + (len(e.rtt_sequence) if e.rtt_sequence else 0)

                gb = []
                gb.append(f"LOCUS       {locus_name:<16s} {seq_len} bp    DNA     linear   SYN")
                gb.append(f"DEFINITION  pegRNA targeting {gene}, {e.edit_type or 'unknown'} edit.")
                if e.paper and e.paper.pmid:
                    gb.append(f"ACCESSION   PMID:{e.paper.pmid}")
                gb.append(f"COMMENT     Editing efficiency: {e.editing_efficiency or 'N/A'}%")
                if e.prime_editor:
                    gb.append(f"COMMENT     Prime editor: {e.prime_editor}")
                if e.cell_type:
                    gb.append(f"COMMENT     Cell type: {e.cell_type}")
                gb.append("FEATURES             Location/Qualifiers")
                gb.append(f'     misc_feature    1..{spacer_end}')
                gb.append(f'                     /label="Spacer"')
                gb.append(f'                     /note="Spacer sequence targeting {gene}"')
                gb.append(f'     misc_feature    {spacer_end+1}..{scaffold_end}')
                gb.append(f'                     /label="Scaffold"')
                gb.append(f'                     /note="SpCas9 tracrRNA scaffold"')
                if e.rtt_sequence:
                    gb.append(f'     misc_feature    {scaffold_end+1}..{rtt_end}')
                    gb.append(f'                     /label="RTT"')
                    gb.append(f'                     /note="Reverse transcriptase template ({len(e.rtt_sequence)}nt)"')
                if e.pbs_sequence:
                    pbs_start = rtt_end + 1 if e.rtt_sequence else scaffold_end + 1
                    pbs_end = pbs_start + len(e.pbs_sequence) - 1
                    gb.append(f'     misc_feature    {pbs_start}..{pbs_end}')
                    gb.append(f'                     /label="PBS"')
                    gb.append(f'                     /note="Primer binding site ({len(e.pbs_sequence)}nt)"')
                gb.append("ORIGIN")
                # Format sequence in GenBank style (60 chars per line, numbered)
                for i in range(0, len(full), 60):
                    chunk = full[i:i+60].lower()
                    parts = " ".join([chunk[j:j+10] for j in range(0, len(chunk), 10)])
                    gb.append(f"{i+1:>9} {parts}")
                gb.append("//")
                gb_records.append("\n".join(gb))

            if gb_records:
                st.download_button(
                    f"Download GenBank/SnapGene ({len(gb_records)} seqs)",
                    "\n".join(gb_records),
                    file_name="pegrna_search_results.gb",
                    mime="text/plain",
                )
            else:
                st.caption("No full sequences available for GenBank export")

        # Detail expander for selected entry
        selected_id = st.selectbox("View entry details", [r["ID"] for r in rows])
        if selected_id:
            entry = session.query(PegRNAEntry).get(selected_id)
            if entry:
                with st.expander("Entry Details", expanded=True):
                    # Full pegRNA construct display (most prominent)
        
                    # Build extension from components if not stored
                    ext_seq = getattr(entry, 'extension_sequence', None)
                    if not ext_seq and entry.rtt_sequence and entry.pbs_sequence:
                        ext_seq = entry.rtt_sequence + entry.pbs_sequence

                    if entry.full_sequence:
                        st.markdown("**Full pegRNA Sequence**")
                        st.code(entry.full_sequence, language=None)
                        full = entry.full_sequence
                        idx = full.find(SCAFFOLD)
                        if idx >= 0:
                            spacer_part = full[:idx]
                            ext_part = full[idx + len(SCAFFOLD):]
                            st.caption(
                                f"Structure: 5'-[Spacer: {len(spacer_part)}nt] "
                                f"[Scaffold: {len(SCAFFOLD)}nt] "
                                f"[3' Extension: {len(ext_part)}nt]-3'"
                            )
                    elif entry.spacer_sequence and ext_seq:
                        reconstructed = entry.spacer_sequence + SCAFFOLD + ext_seq
                        st.markdown("**Full pegRNA Sequence** *(reconstructed)*")
                        st.code(reconstructed, language=None)
                        st.caption(
                            f"Structure: 5'-[Spacer: {len(entry.spacer_sequence)}nt] "
                            f"[Scaffold: {len(SCAFFOLD)}nt] "
                            f"[3' Extension: {len(ext_seq)}nt]-3'"
                        )

                    # 3' Extension (repair template)
                    if ext_seq:
                        st.markdown("**3' Extension (RTT + PBS)**")
                        st.code(ext_seq, language=None)
                        rtt_len = entry.rtt_length or (len(entry.rtt_sequence) if entry.rtt_sequence else "?")
                        pbs_len = entry.pbs_length or (len(entry.pbs_sequence) if entry.pbs_sequence else "?")
                        st.caption(f"RTT: {rtt_len}nt + PBS: {pbs_len}nt = {len(ext_seq)}nt total")

                    dcol1, dcol2 = st.columns(2)
                    with dcol1:
                        st.markdown("**Individual Sequences**")
                        st.code(f"Spacer:  {entry.spacer_sequence or 'N/A'}")
                        st.code(f"PBS:     {entry.pbs_sequence or 'N/A'} ({entry.pbs_length or '?'}nt)")
                        st.code(f"RTT:     {entry.rtt_sequence or 'N/A'} ({entry.rtt_length or '?'}nt)")
                        if entry.three_prime_extension:
                            st.caption(f"3' structural motif: {entry.three_prime_extension}")
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

                    # ClinVar matches for this entry
                    if _has_clinvar and hasattr(entry, 'clinvar_matches') and entry.clinvar_matches:
                        st.markdown("**ClinVar Matches**")
                        cv_rows = []
                        for m in entry.clinvar_matches[:10]:
                            cv = m.clinvar_variant
                            sig = cv.clinical_significance or "Unknown"
                            cv_rows.append({
                                "Variant": cv.name or f"{cv.gene_symbol} {cv.hgvs_cdna or ''}",
                                "Significance": sig,
                                "Condition": (cv.phenotype_list or "N/A")[:80],
                                "Confidence": f"{m.match_confidence:.2f}",
                                "Match Type": m.match_type,
                                "Link": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{cv.variation_id}/"
                            })
                        cv_df = pd.DataFrame(cv_rows)

                        # Color-code significance
                        def color_sig(val):
                            if "pathogenic" in val.lower() and "benign" not in val.lower():
                                return "color: #ff4444"
                            elif "benign" in val.lower():
                                return "color: #44aa44"
                            return "color: #888888"

                        styled = cv_df.style.map(color_sig, subset=["Significance"])
                        st.dataframe(styled, use_container_width=True)

                    if entry.paper:
                        st.markdown("**Source**")
                        st.write(strip_html(entry.paper.title or ""))
                        st.write(strip_html(entry.paper.authors or ""))
                        st.write(f"{entry.paper.journal or 'N/A'} ({entry.paper.year or 'N/A'})")
                        links = []
                        if entry.paper.pmid:
                            links.append(f"[PubMed](https://pubmed.ncbi.nlm.nih.gov/{entry.paper.pmid}/)")
                        if entry.paper.pmcid:
                            links.append(f"[PMC](https://pmc.ncbi.nlm.nih.gov/articles/{entry.paper.pmcid}/)")
                        if entry.paper.doi:
                            links.append(f"[DOI](https://doi.org/{entry.paper.doi})")
                        if links:
                            st.markdown(" | ".join(links))
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
    col3.metric("Total Entries", f"{stats['total_entries']:,}")
    col4.metric("Validated Entries", stats["validated_entries"])

    col5, col6, col7 = st.columns(3)
    col5.metric("Unique Genes", stats["unique_genes"])
    col6.metric("Unique Cell Types", stats["unique_cell_types"])
    col7.metric("Unique Organisms", stats["unique_organisms"])

    # Charts - using SQL aggregations for performance (600k+ entries)
    if stats["total_entries"] > 0:
        from sqlalchemy import func, distinct, case

        # --- Edit type distribution ---
        st.subheader("Distribution of Edit Types")
        VALID_EDIT_TYPES = {
            "substitution", "insertion", "deletion", "1bpreplacement",
            "multibpreplacement", "mnv", "onv", "duplication",
            "inversion", "translocation",
        }
        edit_type_counts = (
            session.query(PegRNAEntry.edit_type, func.count(PegRNAEntry.id))
            .filter(PegRNAEntry.edit_type.isnot(None))
            .group_by(PegRNAEntry.edit_type)
            .all()
        )
        if edit_type_counts:
            et_data = []
            other_count = 0
            for et, cnt in edit_type_counts:
                if et.lower().strip() in VALID_EDIT_TYPES:
                    et_data.append({"edit_type": et, "count": cnt})
                else:
                    other_count += cnt
            if other_count > 0:
                et_data.append({"edit_type": "Other", "count": other_count})
            et_df = pd.DataFrame(et_data).sort_values("count", ascending=False)
            fig = px.pie(
                et_df, names="edit_type", values="count", title="Edit Types",
                hole=0.3,
            )
            fig.update_traces(textposition="inside", textinfo="percent+label")
            fig.update_layout(showlegend=False)
            st.plotly_chart(fig, use_container_width=True)

        # --- Efficiency distribution ---
        st.subheader("Editing Efficiency Distribution")
        eff_rows = (
            session.query(PegRNAEntry.editing_efficiency)
            .filter(
                PegRNAEntry.editing_efficiency.isnot(None),
                PegRNAEntry.editing_efficiency >= 0,
                PegRNAEntry.editing_efficiency <= 100,
            )
            .limit(100000)
            .all()
        )
        if eff_rows:
            eff_df = pd.DataFrame(eff_rows, columns=["efficiency"])
            fig = px.histogram(
                eff_df, x="efficiency", nbins=50,
                title=f"Editing Efficiency Distribution (%) - {len(eff_rows):,} entries",
            )
            fig.update_layout(xaxis_title="Editing Efficiency (%)", yaxis_title="Count")
            st.plotly_chart(fig, use_container_width=True)

        # --- Top genes ---
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

        # --- Prime editor usage ---
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

        # --- pegRNA vs epegRNA efficiency ---
        st.subheader("pegRNA vs epegRNA Efficiency")
        type_eff_rows = (
            session.query(PegRNAEntry.editing_efficiency, PegRNAEntry.pegrna_type)
            .filter(
                PegRNAEntry.editing_efficiency.isnot(None),
                PegRNAEntry.pegrna_type.in_(["pegRNA", "epegRNA"]),
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
            fig.update_layout(xaxis_title="pegRNA Type", yaxis_title="Editing Efficiency (%)")
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No entries with pegRNA/epegRNA type classification and efficiency data.")

        # --- Papers by year ---
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

        # ====================
        # NEW CHARTS
        # ====================

        # --- Delivery Methods ---
        st.subheader("Delivery Methods")
        delivery_counts = (
            session.query(PegRNAEntry.delivery_method, func.count(PegRNAEntry.id).label("count"))
            .filter(PegRNAEntry.delivery_method.isnot(None))
            .group_by(PegRNAEntry.delivery_method)
            .order_by(func.count(PegRNAEntry.id).desc())
            .limit(15)
            .all()
        )
        if delivery_counts:
            del_df = pd.DataFrame(delivery_counts, columns=["Method", "Count"])
            fig = px.bar(del_df, x="Method", y="Count", title="Top 15 Delivery Methods")
            st.plotly_chart(fig, use_container_width=True)

        # --- Organism Distribution ---
        st.subheader("Organism Distribution")
        org_counts = (
            session.query(PegRNAEntry.target_organism, func.count(PegRNAEntry.id).label("count"))
            .filter(PegRNAEntry.target_organism.isnot(None))
            .group_by(PegRNAEntry.target_organism)
            .order_by(func.count(PegRNAEntry.id).desc())
            .limit(15)
            .all()
        )
        if org_counts:
            org_df = pd.DataFrame(org_counts, columns=["Organism", "Count"])
            fig = px.bar(org_df, x="Organism", y="Count", title="Top 15 Target Organisms")
            st.plotly_chart(fig, use_container_width=True)

        # --- PBS Length Distribution (using SQL aggregation) ---
        st.subheader("PBS Length Distribution")
        pbs_agg = (
            session.query(PegRNAEntry.pbs_length, func.count(PegRNAEntry.id).label("count"))
            .filter(PegRNAEntry.pbs_length.isnot(None), PegRNAEntry.pbs_length > 0, PegRNAEntry.pbs_length <= 50)
            .group_by(PegRNAEntry.pbs_length)
            .order_by(PegRNAEntry.pbs_length)
            .all()
        )
        if pbs_agg:
            pbs_df = pd.DataFrame(pbs_agg, columns=["PBS Length (nt)", "Count"])
            total_pbs = pbs_df["Count"].sum()
            fig = px.bar(
                pbs_df, x="PBS Length (nt)", y="Count",
                title=f"PBS Length Distribution ({total_pbs:,} entries)",
            )
            st.plotly_chart(fig, use_container_width=True)

        # --- RTT Length Distribution (using SQL aggregation) ---
        st.subheader("RTT Length Distribution")
        rtt_agg = (
            session.query(PegRNAEntry.rtt_length, func.count(PegRNAEntry.id).label("count"))
            .filter(PegRNAEntry.rtt_length.isnot(None), PegRNAEntry.rtt_length > 0, PegRNAEntry.rtt_length <= 100)
            .group_by(PegRNAEntry.rtt_length)
            .order_by(PegRNAEntry.rtt_length)
            .all()
        )
        if rtt_agg:
            rtt_df = pd.DataFrame(rtt_agg, columns=["RTT Length (nt)", "Count"])
            total_rtt = rtt_df["Count"].sum()
            fig = px.bar(
                rtt_df, x="RTT Length (nt)", y="Count",
                title=f"RTT Length Distribution ({total_rtt:,} entries)",
            )
            st.plotly_chart(fig, use_container_width=True)

        # --- Data Completeness Dashboard ---
        st.subheader("Data Completeness")
        total = stats["total_entries"]
        fields = [
            ("spacer_sequence", "Spacer"),
            ("pbs_sequence", "PBS Seq"),
            ("pbs_length", "PBS Length"),
            ("rtt_sequence", "RTT Seq"),
            ("rtt_length", "RTT Length"),
            ("extension_sequence", "3' Extension"),
            ("full_sequence", "Full Seq"),
            ("target_gene", "Gene"),
            ("target_organism", "Organism"),
            ("edit_type", "Edit Type"),
            ("prime_editor", "PE Version"),
            ("cell_type", "Cell Type"),
            ("delivery_method", "Delivery"),
            ("editing_efficiency", "Efficiency"),
            ("nicking_sgrna_seq", "Nick sgRNA"),
            ("three_prime_extension", "3' Extension"),
        ]
        completeness = []
        for col_name, label in fields:
            col_attr = getattr(PegRNAEntry, col_name)
            non_null = session.query(func.count(PegRNAEntry.id)).filter(col_attr.isnot(None)).scalar()
            pct = (non_null / total * 100) if total > 0 else 0
            completeness.append({"Field": label, "Completeness (%)": round(pct, 1)})

        comp_df = pd.DataFrame(completeness)
        # Color: red (<25%), orange (<50%), yellow (<75%), green (>=75%)
        colors = []
        for pct in comp_df["Completeness (%)"]:
            if pct < 25:
                colors.append("#ff4444")
            elif pct < 50:
                colors.append("#ff8800")
            elif pct < 75:
                colors.append("#ffcc00")
            else:
                colors.append("#44aa44")

        fig = go.Figure(go.Bar(
            x=comp_df["Field"], y=comp_df["Completeness (%)"],
            marker_color=colors,
            text=[f"{v}%" for v in comp_df["Completeness (%)"]],
            textposition="outside",
        ))
        fig.update_layout(title="Data Field Completeness (%)", yaxis_range=[0, 105])
        st.plotly_chart(fig, use_container_width=True)

        # --- Gene × Cell Type Heatmap ---
        st.subheader("Gene × Cell Type Heatmap")
        # Get top 10 genes and top 10 cell types
        top_genes = (
            session.query(PegRNAEntry.target_gene)
            .filter(PegRNAEntry.target_gene.isnot(None))
            .group_by(PegRNAEntry.target_gene)
            .order_by(func.count(PegRNAEntry.id).desc())
            .limit(10)
            .all()
        )
        top_cells = (
            session.query(PegRNAEntry.cell_type)
            .filter(PegRNAEntry.cell_type.isnot(None))
            .group_by(PegRNAEntry.cell_type)
            .order_by(func.count(PegRNAEntry.id).desc())
            .limit(10)
            .all()
        )
        if top_genes and top_cells:
            gene_names = [g[0] for g in top_genes]
            cell_names = [c[0] for c in top_cells]

            heatmap_data = (
                session.query(
                    PegRNAEntry.target_gene,
                    PegRNAEntry.cell_type,
                    func.count(PegRNAEntry.id),
                )
                .filter(
                    PegRNAEntry.target_gene.in_(gene_names),
                    PegRNAEntry.cell_type.in_(cell_names),
                )
                .group_by(PegRNAEntry.target_gene, PegRNAEntry.cell_type)
                .all()
            )

            # Build matrix
            matrix = pd.DataFrame(0, index=gene_names, columns=cell_names)
            for gene, cell, count in heatmap_data:
                if gene in matrix.index and cell in matrix.columns:
                    matrix.loc[gene, cell] = count

            fig = px.imshow(
                matrix.values,
                x=matrix.columns.tolist(),
                y=matrix.index.tolist(),
                title="Entry Count: Top 10 Genes × Top 10 Cell Types",
                labels=dict(x="Cell Type", y="Gene", color="Entries"),
                color_continuous_scale="YlOrRd",
                aspect="auto",
            )
            st.plotly_chart(fig, use_container_width=True)

        # ====================
        # ClinVar Statistics (only if data loaded)
        # ====================
        if has_clinvar_data(session):
            st.markdown("---")
            st.header("ClinVar Clinical Variant Data")

            from database.models import ClinVarVariant, PegRNAClinVarMatch

            cv_total = session.query(ClinVarVariant).count()
            matched_entries = session.query(func.count(distinct(PegRNAClinVarMatch.pegrna_entry_id))).scalar()
            matched_variants = session.query(func.count(distinct(PegRNAClinVarMatch.clinvar_variant_id))).scalar()

            cv1, cv2, cv3 = st.columns(3)
            cv1.metric("ClinVar Variants", f"{cv_total:,}")
            cv2.metric("Entries with ClinVar Match", f"{matched_entries:,}")
            cv3.metric("Matched Variants", f"{matched_variants:,}")

            # Clinical significance breakdown
            st.subheader("Clinical Significance of Matched Variants")
            sig_counts = (
                session.query(ClinVarVariant.clinical_significance, func.count(ClinVarVariant.id))
                .join(PegRNAClinVarMatch, PegRNAClinVarMatch.clinvar_variant_id == ClinVarVariant.id)
                .group_by(ClinVarVariant.clinical_significance)
                .order_by(func.count(ClinVarVariant.id).desc())
                .limit(10)
                .all()
            )
            if sig_counts:
                sig_df = pd.DataFrame(sig_counts, columns=["Significance", "Count"])
                fig = px.bar(sig_df, x="Significance", y="Count",
                             title="Clinical Significance of Matched ClinVar Variants")
                st.plotly_chart(fig, use_container_width=True)

            # Top genes with pathogenic matches
            st.subheader("Top Genes with Pathogenic Variant Matches")
            patho_gene_counts = (
                session.query(PegRNAEntry.target_gene, func.count(distinct(PegRNAClinVarMatch.clinvar_variant_id)))
                .join(PegRNAClinVarMatch, PegRNAClinVarMatch.pegrna_entry_id == PegRNAEntry.id)
                .join(ClinVarVariant, PegRNAClinVarMatch.clinvar_variant_id == ClinVarVariant.id)
                .filter(ClinVarVariant.clinical_significance.ilike("%pathogenic%"))
                .group_by(PegRNAEntry.target_gene)
                .order_by(func.count(distinct(PegRNAClinVarMatch.clinvar_variant_id)).desc())
                .limit(15)
                .all()
            )
            if patho_gene_counts:
                pg_df = pd.DataFrame(patho_gene_counts, columns=["Gene", "Pathogenic Variants"])
                fig = px.bar(pg_df, x="Gene", y="Pathogenic Variants",
                             title="Genes with Most Pathogenic ClinVar Variant Matches")
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
        "Enter PMIDs, DOIs, or **any publisher URL** to extract pegRNA data. "
        "One identifier per line. Supports Nature, MDPI, Cell, Wiley, Science, "
        "bioRxiv, Springer, Frontiers, PLOS, and more."
    )

    input_text = st.text_area(
        "Paper Identifiers",
        height=150,
        placeholder="34608327\nhttps://www.mdpi.com/2073-4425/13/12/2348\nhttps://www.nature.com/articles/s41586-024-07259-6\n10.1038/s41587-022-01613-7",
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
                    elif parsed["type"] == "url":
                        log_area.info(f"Fetching DOI from URL: {ident[:80]}...")
                        pmid = resolve_to_pmid(parsed)
                        if pmid:
                            pmids.append(pmid)
                        else:
                            log_area.warning(f"Could not extract DOI from URL: {ident}")
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

    st.markdown(
        "Export the full database or use **Search & Browse** to filter and export a subset. "
        "Search results can be exported as CSV, FASTA, or SnapGene-compatible GenBank format."
    )

    export_format = st.selectbox("Format", ["CSV", "JSON"])

    from sqlalchemy import func as _func
    entry_count = session.query(_func.count(PegRNAEntry.id)).scalar()

    st.write(f"**{entry_count:,} entries** to export")

    if entry_count > 0 and st.button("Generate Export"):
        with st.spinner(f"Generating {export_format} export for {entry_count:,} entries..."):
            from sqlalchemy import text
            engine = session.get_bind()

            sql = """
                SELECT
                    e.id, e.entry_name, e.pegrna_type,
                    e.spacer_sequence, e.pbs_sequence, e.pbs_length,
                    e.rtt_sequence, e.rtt_length, e.extension_sequence,
                    e.three_prime_extension, e.full_sequence, e.nicking_sgrna_seq,
                    e.target_gene, e.target_locus, e.target_organism,
                    e.edit_type, e.edit_description, e.intended_mutation,
                    e.prime_editor, e.cell_type, e.delivery_method,
                    e.editing_efficiency, e.product_purity, e.indel_frequency,
                    e.confidence_score, e.validated,
                    p.pmid AS paper_pmid, p.doi AS paper_doi, p.year AS paper_year
                FROM pegrna_entries e
                LEFT JOIN papers p ON e.paper_id = p.id
            """
            df = pd.read_sql(text(sql), engine)

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

## ClinVar Integration

Entries can be cross-referenced with NCBI ClinVar to identify pegRNAs that target
clinically relevant genetic variants. Matches are scored by gene, edit type alignment,
and mutation description overlap.

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
