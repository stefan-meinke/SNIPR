from Bio import SeqIO

def load_genome(fasta_path):
    """
    Load genome as dictionary of Seq objects
    """
    return SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
