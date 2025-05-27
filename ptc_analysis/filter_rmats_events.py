import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Path to raw rMATS output (*.MATS.JC.txt)')
    parser.add_argument('--output', required=True, help='Path to save filtered file')
    parser.add_argument('--fdr', type=float, default=0.01, help='FDR cutoff')
    parser.add_argument('--dpsi', type=float, default=0.15, help='Minimum absolute IncLevelDifference')
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t')

    filtered = df[
        (df['FDR'] < args.fdr) &
        (df['IncLevelDifference'].abs() >= args.dpsi)
    ]

    filtered.to_csv(args.output, sep='\t', index=False)
    print(f"✅ Filtered {len(filtered)} events saved to: {args.output}")

if __name__ == "__main__":
    main()
