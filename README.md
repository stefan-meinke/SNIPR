
# SNIPR (Splicing-Related NMD Impact and Protein Regulation)

A modular Python pipeline to analyze the impact of alternative splicing on coding sequence integrity and gene expression. It simulates spliced isoforms, translates coding sequences, detects premature stop codons (PTCs), and predicts likelihood of nonsense-mediated decay (NMD).

---

## Features

- Processes [rMATS-turbo](https://github.com/Xinglab/rmats-turbo) ([Wang et al., 2024](https://www.nature.com/articles/s41596-023-00944-2)) output (SE, RI, MXE, A3SS, A5SS)
- Reconstructs spliced transcripts using genome FASTA + GTF
- Simulates exon skipping or inclusion
- Detects ORF truncation and premature stop codons (PTCs)
- Predicts NMD-triggering events using the 50-nt rule
- Outputs detailed summaries per transcript and splicing event
- Supports batch analysis across multiple datasets

---

## Running via Guix

### Step 1: Launch a Guix shell

```bash
guix shell -m manifest.scm
```

### Step 2: Install the tool inside the Guix shell

```bash
git clone https://github.com/stefan-meinke/SNIPR.git
cd snipr
pip install .
```

Now you can use the CLI command `ptc_analysis`.

---

## Inputs

- `*.MATS.JC.txt` and `fromGTF.*.txt` from rMATS
- Reference GTF (e.g., Ensembl)
- Genome FASTA file

### Example Structure

```
datasets/
└── dataset1/
    ├── SE.MATS.JC.txt
    ├── fromGTF.SE.txt
    ├── RI.MATS.JC.txt
    └── fromGTF.RI.txt
```

---

## Filtering rMATS Events

Before analysis, filter events for significance (defaults: fdr < 0.01, ):

```bash
python -m ptc_analysis.filter_rmats_events \
  --input datasets/dataset1/SE.MATS.JC.txt \
  --output datasets/dataset1/SE.MATS.JC.filtered.txt \
  --fdr 0.01 \
  --dpsi 0.15
```

---

## Usage

### Analyze a single splice type (e.g. skipped exon (SE) events:

splice types are:
- SE (skipped exons)
- RI (retained introns)
- MXE (mutually exclusive exons)
- A3SS (alternative 3' splice site)
- A5SS (alternative 5' splice site)

```bash
ptc_analysis \
  --dataset_dir datasets/dataset1 \
  --splice_type SE \ 
  --output_dir results/dataset1/SE \
  --gtf reference.gtf \
  --fasta genome.fa
```

### Run on all splice types for one dataset:

```bash
python run_ptc_batch.py \
  --single_dataset datasets/dataset1 \
  --output_dir results \
  --gtf reference.gtf \
  --fasta genome.fa
```

### Run on multiple datasets:

```bash
python run_ptc_batch.py \
  --dataset_dir datasets \
  --output_dir results \
  --gtf reference.gtf \
  --fasta genome.fa
```

## Extract Sequences of Spliced Exons and 250 bp upstream and downstream

```bash
python extract_flanks.py \
  --dataset datasets/dataset/SE.MATS.JC.filtered.txt \
  --fromgtf datasets/dataset/fromGTF.SE.ftxt \
  --genome genome.fa \
  --splice_type SE \
  --output spliced_exons_with_flanks.fa
```

---

## Output

- `orf_disruption_results.csv`: per-transcript disruption and NMD prediction
- `skipped_transcripts_log.csv`: transcripts where simulation failed
- `spliced_exons_with_flanks.fa`:
FASTA file of flanked spliced regions

---

## Requirements (handled by Guix)

```scheme
(specifications->manifest
  (list
    "python"
    "python-pandas"
    "python-biopython"
    "python-gffutils"
    "python-setuptools"
    "python-pip"
    "coreutils"
    "grep"
    "bash"
  ))
```

---

## License

MIT License

---

## Contact

Author: Stefan Meinke 
GitHub: [stefan-meinke/snipr](https://github.com/stefan-meinke/SNIPR.git)

---

## Citation

> SNIPR: A tool to predict coding disruption and NMD from alternative splicing. GitHub. 2025.
