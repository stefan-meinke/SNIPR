import os
import argparse
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Run ptc_analysis on all splice types")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--dataset_dir', help='Directory containing multiple datasets (each subfolder is a dataset)')
    group.add_argument('--single_dataset', help='Path to a single dataset directory')

    parser.add_argument('--output_dir', required=True, help='Base directory to store results')
    parser.add_argument('--gtf', required=True, help='Reference GTF file')
    parser.add_argument('--fasta', required=True, help='Genome FASTA file')

    args = parser.parse_args()
    splice_types = ["SE", "RI", "MXE", "A3SS", "A5SS"]

    # Determine datasets to process
    if args.dataset_dir:
        datasets = [(ds, os.path.join(args.dataset_dir, ds))
                    for ds in os.listdir(args.dataset_dir)
                    if os.path.isdir(os.path.join(args.dataset_dir, ds))]
    else:
        ds_name = os.path.basename(os.path.normpath(args.single_dataset))
        datasets = [(ds_name, args.single_dataset)]

    for ds_name, dataset_path in datasets:
        for sp in splice_types:
            result_dir = os.path.join(args.output_dir, ds_name, sp)
            result_file = os.path.join(result_dir, "orf_disruption_results.csv")

            if os.path.exists(result_file):
                print(f"✅ Skipping {ds_name} | {sp} (already processed)")
                continue

            print(f"🔄 Processing {ds_name} | {sp}")
            os.makedirs(result_dir, exist_ok=True)

            subprocess.run([
                "ptc_analysis",
                "--dataset_dir", dataset_path,
                "--splice_type", sp,
                "--output_dir", result_dir,
                "--gtf", args.gtf,
                "--fasta", args.fasta
            ])

if __name__ == "__main__":
    main()
