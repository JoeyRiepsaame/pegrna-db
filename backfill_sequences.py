"""Back-fill extension_sequence and full_sequence for existing entries.

Run: python backfill_sequences.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import config
from database.models import init_db, PegRNAEntry
from sqlalchemy import text

# Common SpCas9 scaffold (most used in prime editing)
SCAFFOLD = "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"


def main():
    Session = init_db(str(config.DATABASE_PATH))
    session = Session()

    # Step 0: Ensure the extension_sequence column exists
    try:
        session.execute(text("SELECT extension_sequence FROM pegrna_entries LIMIT 1"))
        print("extension_sequence column exists")
    except Exception:
        print("Adding extension_sequence column...")
        session.execute(text("ALTER TABLE pegrna_entries ADD COLUMN extension_sequence TEXT"))
        session.commit()

    # Step 1: Back-fill extension_sequence = RTT + PBS for entries that have both
    print("Back-filling extension_sequence from RTT + PBS...")
    result = session.execute(text("""
        UPDATE pegrna_entries
        SET extension_sequence = rtt_sequence || pbs_sequence
        WHERE rtt_sequence IS NOT NULL
          AND pbs_sequence IS NOT NULL
          AND (extension_sequence IS NULL OR extension_sequence = '')
    """))
    ext_updated = result.rowcount
    session.commit()
    print(f"  Updated {ext_updated} entries with extension_sequence")

    # Step 2: Reconstruct full_sequence = spacer + scaffold + RTT + PBS
    print("Reconstructing full_sequence from spacer + scaffold + extension...")
    scaffold = SCAFFOLD
    result = session.execute(text(f"""
        UPDATE pegrna_entries
        SET full_sequence = spacer_sequence || '{scaffold}' || rtt_sequence || pbs_sequence
        WHERE spacer_sequence IS NOT NULL
          AND rtt_sequence IS NOT NULL
          AND pbs_sequence IS NOT NULL
          AND (full_sequence IS NULL OR full_sequence = '')
    """))
    full_updated = result.rowcount
    session.commit()
    print(f"  Updated {full_updated} entries with full_sequence")

    # Step 3: Decompose existing full_sequences into components
    print("Decomposing existing full_sequences into components...")
    entries_with_full = (
        session.query(PegRNAEntry)
        .filter(PegRNAEntry.full_sequence.isnot(None))
        .filter(
            (PegRNAEntry.spacer_sequence.is_(None)) |
            (PegRNAEntry.extension_sequence.is_(None))
        )
        .all()
    )
    decomposed = 0
    for entry in entries_with_full:
        full = entry.full_sequence
        for sc in [
            SCAFFOLD,
            "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
            "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTG",
        ]:
            idx = full.find(sc)
            if idx >= 0:
                if 15 <= idx <= 25:
                    if not entry.spacer_sequence:
                        entry.spacer_sequence = full[:idx]
                    ext = full[idx + len(sc):]
                    if ext and not entry.extension_sequence:
                        entry.extension_sequence = ext
                    decomposed += 1
                break
    if decomposed:
        session.commit()
    print(f"  Decomposed {decomposed} full_sequences")

    # Print final statistics
    total = session.execute(text("SELECT COUNT(*) FROM pegrna_entries")).scalar()
    ext_count = session.execute(text(
        "SELECT COUNT(*) FROM pegrna_entries WHERE extension_sequence IS NOT NULL"
    )).scalar()
    full_count = session.execute(text(
        "SELECT COUNT(*) FROM pegrna_entries WHERE full_sequence IS NOT NULL"
    )).scalar()
    spacer_count = session.execute(text(
        "SELECT COUNT(*) FROM pegrna_entries WHERE spacer_sequence IS NOT NULL"
    )).scalar()

    print(f"\nFinal statistics:")
    print(f"  Total entries: {total:,}")
    print(f"  Spacer:    {spacer_count:,} ({100*spacer_count/total:.1f}%)")
    print(f"  Extension: {ext_count:,} ({100*ext_count/total:.1f}%)")
    print(f"  Full seq:  {full_count:,} ({100*full_count/total:.1f}%)")

    session.close()


if __name__ == "__main__":
    main()
