"""Fix missing efficiency data for PMC12008803 (PMID 40120586).

Paper: "High-throughput screening of human genetic variants by pooled prime editing"
21,290 entries across 4 supplementary files, all missing editing_efficiency.

Processing order (verified by DB entry ID order):
  mmc4 (10,954 rows) - SMARCB1 screen → D10 R1/R2 average
  mmc6 (2,696 rows)  - MLH1 exon 10 screen → PEmax_D20
  mmc8 (3,748 rows)  - MLH1 non-coding screen → PEmax_D20
  mmc2 (3,892 rows)  - Scaffold comparison → percentage_editing
"""
import numpy as np
import pandas as pd
from database.models import init_db, PegRNAEntry

BASE = "data/raw_papers/PMC12008803/PMC12008803"
PAPER_ID = 44


def load_efficiencies():
    """Load efficiency values from each Excel file in DB entry order."""
    results = []

    # 1. mmc4 - SMARCB1 (10,954 rows)
    # Use average of D10_R1 and D10_R2 (Day 10, standard PE screening timepoint)
    df = pd.read_excel(f"{BASE}/mmc4.xlsx", header=1)
    print(f"mmc4: {len(df)} rows")
    for _, row in df.iterrows():
        r1 = row.get("D10_R1_percentage_editing")
        r2 = row.get("D10_R2_percentage_editing")
        vals = [v for v in [r1, r2] if pd.notna(v)]
        eff = np.mean(vals) if vals else None
        results.append({
            "efficiency": round(eff, 2) if eff is not None else None,
            "prime_editor": "PEmax",
            "cell_type": "HAP1",
            "source": "mmc4",
        })

    # 2. mmc6 - MLH1 exon 10 (2,696 rows)
    df = pd.read_excel(f"{BASE}/mmc6.xlsx", header=1)
    print(f"mmc6: {len(df)} rows")
    for _, row in df.iterrows():
        eff = row.get("PEmax_D20_percentage_editing")
        results.append({
            "efficiency": round(float(eff), 2) if pd.notna(eff) else None,
            "prime_editor": "PEmax",
            "cell_type": "HAP1",
            "source": "mmc6",
        })

    # 3. mmc8 - MLH1 non-coding (3,748 rows)
    df = pd.read_excel(f"{BASE}/mmc8.xlsx", header=1)
    print(f"mmc8: {len(df)} rows")
    for _, row in df.iterrows():
        eff = row.get("PEmax_D20_percentage_editing")
        results.append({
            "efficiency": round(float(eff), 2) if pd.notna(eff) else None,
            "prime_editor": "PEmax",
            "cell_type": "HAP1",
            "source": "mmc8",
        })

    # 4. mmc2 - Scaffold comparison (3,892 rows)
    df = pd.read_excel(f"{BASE}/mmc2.xlsx", header=1)
    print(f"mmc2: {len(df)} rows")
    for _, row in df.iterrows():
        eff = row.get("percentage_editing")
        results.append({
            "efficiency": round(float(eff), 2) if pd.notna(eff) else None,
            "prime_editor": "PEmax",
            "cell_type": "HEK293T",
            "source": "mmc2",
        })

    return results


def verify_alignment(session, entries, excel_data):
    """Verify that Excel rows align with DB entries by spot-checking spacer sequences."""
    checks = [0, 100, 1000, 10953, 10954, 13649, 17397, 21289]  # boundary indices
    all_ok = True
    for idx in checks:
        if idx >= len(entries) or idx >= len(excel_data):
            continue
        db_spacer = entries[idx].spacer_sequence
        # Load the spacer from the appropriate file
        source = excel_data[idx]["source"]
        # Offset within each file
        offsets = {"mmc4": 0, "mmc6": 10954, "mmc8": 10954 + 2696, "mmc2": 10954 + 2696 + 3748}
        file_idx = idx - offsets[source]
        fname = f"{BASE}/{source}"
        df = pd.read_excel(fname, header=1, skiprows=range(1, file_idx + 1), nrows=1)
        if len(df) > 0:
            excel_spacer = str(df.iloc[0].get("protospacer", "")).strip()
            match = db_spacer == excel_spacer
            if not match:
                print(f"  MISMATCH at idx {idx}: DB={db_spacer}, Excel={excel_spacer} ({source} row {file_idx})")
                all_ok = False
            else:
                print(f"  OK at idx {idx}: {db_spacer} ({source} row {file_idx})")
    return all_ok


def main():
    print("Loading efficiency data from Excel files...")
    excel_data = load_efficiencies()
    print(f"Total Excel rows: {len(excel_data)}")

    eff_count = sum(1 for d in excel_data if d["efficiency"] is not None)
    print(f"Rows with efficiency values: {eff_count} ({100*eff_count/len(excel_data):.1f}%)")

    Session = init_db("data/pegrna_database.db")
    session = Session()

    # Load all DB entries for this paper, ordered by ID
    entries = (
        session.query(PegRNAEntry)
        .filter_by(paper_id=PAPER_ID)
        .order_by(PegRNAEntry.id)
        .all()
    )
    print(f"DB entries for paper {PAPER_ID}: {len(entries)}")

    if len(entries) != len(excel_data):
        print(f"ERROR: Count mismatch! DB={len(entries)}, Excel={len(excel_data)}")
        session.close()
        return

    # Quick alignment verification using spacer sequences at file boundaries
    print("\nVerifying alignment (spot checks)...")
    import pandas as pd
    files_info = [
        ("mmc4.xlsx", 0, 10954),
        ("mmc6.xlsx", 10954, 2696),
        ("mmc8.xlsx", 10954 + 2696, 3748),
        ("mmc2.xlsx", 10954 + 2696 + 3748, 3892),
    ]
    all_ok = True
    for fname, db_offset, nrows in files_info:
        df = pd.read_excel(f"{BASE}/{fname}", header=1, nrows=3)
        for i in range(min(3, len(df))):
            db_spacer = entries[db_offset + i].spacer_sequence
            xl_spacer = str(df.iloc[i]["protospacer"]).strip()
            ok = "OK" if db_spacer == xl_spacer else "MISMATCH"
            if ok == "MISMATCH":
                all_ok = False
            print(f"  {fname} row {i} (DB idx {db_offset + i}): {ok}")
    if not all_ok:
        print("ALIGNMENT FAILED - aborting!")
        session.close()
        return
    print("  All checks passed.")

    # Update entries
    print("\nUpdating entries...")
    updated_eff = 0
    updated_pe = 0
    updated_ct = 0
    for i, (entry, data) in enumerate(zip(entries, excel_data)):
        if data["efficiency"] is not None and entry.editing_efficiency is None:
            entry.editing_efficiency = data["efficiency"]
            updated_eff += 1
        if data["prime_editor"] and not entry.prime_editor:
            entry.prime_editor = data["prime_editor"]
            updated_pe += 1
        if data["cell_type"] and not entry.cell_type:
            entry.cell_type = data["cell_type"]
            updated_ct += 1
        if not entry.target_organism:
            entry.target_organism = "Human"

    print(f"Updated efficiency: {updated_eff}")
    print(f"Updated prime_editor: {updated_pe}")
    print(f"Updated cell_type: {updated_ct}")

    # Source distribution
    from collections import Counter
    source_counts = Counter()
    eff_by_source = Counter()
    for data in excel_data:
        source_counts[data["source"]] += 1
        if data["efficiency"] is not None:
            eff_by_source[data["source"]] += 1
    print("\nBy source file:")
    for src in ["mmc4", "mmc6", "mmc8", "mmc2"]:
        print(f"  {src}: {source_counts[src]} total, {eff_by_source[src]} with efficiency ({100*eff_by_source[src]/source_counts[src]:.1f}%)")

    session.commit()
    print("\nCommitted to database.")
    session.close()


if __name__ == "__main__":
    main()
