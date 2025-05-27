import os
import pandas as pd

from .data_loader import load_filtered_events
from .transcript_model import build_gtf_db, get_transcripts_for_gene
from .splicing_simulator import extract_transcript_sequence, simulate_exon_skipping
from .orf_analysis import (
    translate_cds,
    check_orf_disruption,
    find_stop_position,
    get_last_cds_junction,
    check_nmd_with_junction
)
from .utils import load_genome


def run_analysis(dataset_dir, splice_type, output_dir, gtf_file, fasta_file):
    """
    Run PTC/NMD analysis for a given rMATS dataset and splice type.
    """
    os.makedirs(output_dir, exist_ok=True)

    filtered_rmats = os.path.join(dataset_dir, f"{splice_type}.MATS.JC.filtered.txt")
    gtf_events = os.path.join(dataset_dir, f"fromGTF.{splice_type}.txt")

    print("Loading genome FASTA...")
    genome = load_genome(fasta_file)

    print("Building or loading GTF database...")
    db = build_gtf_db(gtf_file)

    print(f"Loading filtered {splice_type} splicing events...")
    events = load_filtered_events(filtered_rmats, gtf_events, splice_type)

    results = []
    skipped_transcripts = []

    coord_keys = {
        "SE": ("exonStart_0base", "exonEnd"),
        "RI": ("riExonStart_0base", "riExonEnd"),
        "MXE": ("1stExonStart_0base", "1stExonEnd"),
        "A3SS": ("longExonStart_0base", "longExonEnd"),
        "A5SS": ("longExonStart_0base", "longExonEnd")
    }
    start_key, end_key = coord_keys[splice_type]

    for idx, row in events.iterrows():
        gene_id = row['GeneID']
        exon_start = int(row[start_key]) + 1  # rMATS is 0-based
        exon_end = int(row[end_key])
        skip_coords = (exon_start, exon_end)
        dpsi = float(row['IncLevelDifference'])

        try:
            transcripts = get_transcripts_for_gene(db, gene_id)
        except Exception as e:
            print(f"Skipping gene {gene_id} (error: {e})")
            continue

        for tx in transcripts:
            try:
                ref_seq = extract_transcript_sequence(tx, genome, db)

                try:
                    alt_seq = simulate_exon_skipping(tx, skip_coords, db, genome)
                except ValueError as e:
                    skipped_transcripts.append({
                        "GeneID": gene_id,
                        "TranscriptID": tx.id,
                        "ExonStart": exon_start,
                        "ExonEnd": exon_end,
                        "Reason": str(e)
                    })
                    continue

                ref_prot = translate_cds(ref_seq)
                alt_prot = translate_cds(alt_seq)

                disrupted = check_orf_disruption(ref_prot, alt_prot)

                likely_nmd = False
                if disrupted:
                    ptc_pos = find_stop_position(alt_seq)
                    last_junction_nt = get_last_cds_junction(tx, db)
                    likely_nmd = check_nmd_with_junction(ptc_pos, last_junction_nt)

                result = {
                    "GeneID": gene_id,
                    "TranscriptID": tx.id,
                    "ExonStart": exon_start,
                    "ExonEnd": exon_end,
                    "RefLen_AA": len(ref_prot),
                    "AltLen_AA": len(alt_prot),
                    "TruncationAA": len(ref_prot) - len(alt_prot),
                    "Disrupted": disrupted,
                    "Likely_NMD": likely_nmd,
                    "AS_Direction": "Inclusion" if dpsi > 0 else "Exclusion",
                    "IncLevelDifference": dpsi
                }
                results.append(result)

            except Exception as e:
                print(f"Error processing transcript {tx.id}: {e}")
                continue

    print(f"Finished processing {len(results)} transcript events.")
    results_df = pd.DataFrame(results)
    results_outfile = os.path.join(output_dir, "orf_disruption_results.csv")
    results_df.to_csv(results_outfile, index=False)
    print(f"Results saved to: {results_outfile}")

    if skipped_transcripts:
        skipped_df = pd.DataFrame(skipped_transcripts)
        skipped_file = os.path.join(output_dir, "skipped_transcripts_log.csv")
        skipped_df.to_csv(skipped_file, index=False)
        print(f"Skipped transcripts saved to: {skipped_file}")
    else:
        print("No transcripts skipped due to missing exons.")
