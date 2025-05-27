from Bio.Seq import Seq

def extract_transcript_sequence(transcript, genome, db):
    """
    Reconstruct spliced transcript (inclusion isoform) from CDS exons.
    Returns a nucleotide string.
    """
    exon_seqs = []
    for exon in db.children(transcript, featuretype='CDS', order_by='start'):
        chrom = exon.chrom
        start = exon.start - 1
        end = exon.end
        strand = exon.strand

        raw_seq = genome[chrom][start:end].seq
        seq = raw_seq if strand == '+' else raw_seq.reverse_complement()
        exon_seqs.append(str(seq))

    return ''.join(exon_seqs)


def simulate_exon_skipping(transcript, skip_coords, db, genome):
    """
    Reconstruct transcript with the specified exon skipped.
    skip_coords: (start, end) of exon to skip
    Returns a nucleotide string.
    """
    exon_seqs = []
    found_exon = False

    for exon in db.children(transcript, featuretype='exon', order_by='start'):
        if exon.start == skip_coords[0] and exon.end == skip_coords[1]:
            found_exon = True
            continue  # Skip this exon

        chrom = exon.chrom
        seq = genome[chrom][exon.start - 1:exon.end].seq
        if exon.strand == '-':
            seq = seq.reverse_complement()
        exon_seqs.append(str(seq))

    if not found_exon:
        raise ValueError("Exon to skip not found in transcript")

    return ''.join(exon_seqs)
