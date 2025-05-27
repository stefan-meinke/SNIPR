import pandas as pd

def load_filtered_events(filtered_rmats_file, from_gtf_file, splice_type):
    """
    Load significant splicing events and join with relevant exon coordinates from fromGTF file.
    Coordinates vary depending on splice type (SE, RI, MXE, A3SS, A5SS).
    """

    # Load the filtered rMATS output and fromGTF exon annotations
    filtered = pd.read_csv(filtered_rmats_file, sep='\t')
    gtf_events = pd.read_csv(from_gtf_file, sep='\t')

    # Coordinate columns for each splice type
    coord_map = {
        "SE": ("exonStart_0base", "exonEnd"),
        "RI": ("riExonStart_0base", "riExonEnd"),
        "MXE": ("1stExonStart_0base", "1stExonEnd"),
        "A3SS": ("longExonStart_0base", "longExonEnd"),
        "A5SS": ("longExonStart_0base", "longExonEnd"),
    }

    if splice_type not in coord_map:
        raise ValueError(f"Unsupported splice type: {splice_type}")

    start_col, end_col = coord_map[splice_type]

    # Ensure required columns exist in the GTF input
    missing_cols = [col for col in [start_col, end_col] if col not in gtf_events.columns]
    if missing_cols:
        raise KeyError(f"Missing expected columns in fromGTF file: {missing_cols}")

    # Drop potential duplicates from filtered file before merge
    filtered = filtered.drop(columns=[start_col, end_col], errors='ignore')

    # Merge filtered rMATS results with coordinate info from GTF
    gtf_events = gtf_events[['ID', start_col, end_col]]
    merged = pd.merge(filtered, gtf_events, on='ID', how='inner')

    return merged
