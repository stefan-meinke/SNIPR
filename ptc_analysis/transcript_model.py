import gffutils

def build_gtf_db(gtf_file, db_file="gtf.db"):
    """
    Create or load a gffutils database from the GTF file
    """
    try:
        db = gffutils.FeatureDB(db_file, keep_order=True)
    except:
        db = gffutils.create_db(gtf_file, dbfn=db_file, force=True, keep_order=True, disable_infer_genes=True, disable_infer_transcripts=True)
    return db

def get_transcripts_for_gene(db, gene_id):
    """
    Return all transcript IDs for a gene
    """
    gene = db[gene_id]
    return list(db.children(gene, featuretype='transcript'))
