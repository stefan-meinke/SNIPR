from Bio.Seq import Seq

def translate_cds(mrna_seq):
    """
    Translate mRNA sequence up to the first in-frame stop codon (*).
    """
    return Seq(mrna_seq).translate(to_stop=True)


def check_orf_disruption(ref_protein, alt_protein):
    """
    Returns True if the alternative protein is at least 300 amino acids shorter.
    """
    return len(ref_protein) - len(alt_protein) > 300


def find_stop_position(seq):
    """
    Return nucleotide index (0-based) of first in-frame stop codon.
    If no stop is found, return full length.
    """
    protein = Seq(seq).translate()
    stop_index = protein.find("*")
    if stop_index == -1:
        return len(seq)  # No stop codon; return full CDS length
    return stop_index * 3  # Convert amino acid index to nucleotide


def get_last_cds_junction(transcript, db):
    """
    Get the CDS-relative position (in nucleotides) of the last exon-exon junction in the CDS.
    If only one CDS exon exists, return None (no junctions).
    """
    cds_exons = list(db.children(transcript, featuretype='CDS', order_by='start'))

    if len(cds_exons) < 2:
        return None  # No junctions in CDS

    # Get lengths of each CDS exon
    lengths = [cds.end - cds.start + 1 for cds in cds_exons]

    # Sum up length of all but the last exon to get position of last junction
    last_junction_nt = sum(lengths[:-1])
    return last_junction_nt


def check_nmd_with_junction(ptc_position_nt, last_junction_nt):
    """
    Determines if a premature stop codon (PTC) is >50 nt upstream of the last CDS junction.
    """
    if last_junction_nt is None:
        return False  # No junction = no NMD per canonical rule
    return (last_junction_nt - ptc_position_nt) > 50
