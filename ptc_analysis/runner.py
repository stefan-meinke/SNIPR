import argparse
from .core import run_analysis

def main():
    parser = argparse.ArgumentParser(description="Analyze PTCs and potential NMD from rMATS results")
    parser.add_argument('--dataset_dir', required=True, help='Path to the rMATS output directory')
    parser.add_argument('--splice_type', required=True, choices=['SE', 'RI', 'MXE', 'A3SS', 'A5SS'])
    parser.add_argument('--output_dir', required=True, help='Output directory for results')
    parser.add_argument('--gtf', required=True, help='Path to reference GTF file')
    parser.add_argument('--fasta', required=True, help='Path to genome FASTA file')

    args = parser.parse_args()

    run_analysis(
        dataset_dir=args.dataset_dir,
        splice_type=args.splice_type,
        output_dir=args.output_dir,
        gtf_file=args.gtf,
        fasta_file=args.fasta
    )

if __name__ == "__main__":
    main()