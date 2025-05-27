
# SPLINTR (Splicing-Induced Truncation and Regulation)

A modular Python pipeline to analyze the impact of alternative splicing on coding sequence integrity and gene expression. It simulates spliced isoforms, translates coding sequences, detects premature stop codons (PTCs), and predicts likelihood of nonsense-mediated decay (NMD).

---

## 🚀 Features

- Processes rMATS output (SE, RI, MXE, A3SS, A5SS)
- Reconstructs spliced transcripts using genome FASTA + GTF
- Simulates exon skipping or inclusion
- Detects ORF truncation and premature stop codons (PTCs)
- Predicts NMD-triggering events using the 50-nt rule
- Outputs detailed summaries per transcript and splicing event
- Supports batch analysis across multiple datasets

---

## 🧬 Installation

```bash
git clone https://github.com/yourusername/splintr.git
cd splintr
pip install -e .
```

---

## 📁 Inputs

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

## 🧪 Filtering rMATS Events

Before analysis, filter events for significance:

```bash
python -m ptc_analysis.filter_rmats_events \
  --input SE.MATS.JC.txt \
  --output SE.MATS.JC.filtered.txt \
  --fdr 0.01 \
  --dpsi 0.15
```

---

## 🔧 Usage

### Analyze a single splice type:

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

---

## 📤 Output

- `orf_disruption_results.csv`: per-transcript disruption and NMD prediction
- `skipped_transcripts_log.csv`: transcripts where simulation failed

---

## 📦 Requirements

- Python ≥ 3.12.1
- pandas
- biopython
- gffutils

Install with:

```bash
pip install -r requirements.txt
```

---

## 📄 License

MIT License

---

## 🧠 Contact

Author: Your Name  
GitHub: [yourusername/splintr](https://github.com/yourusername/splintr)

---

## 📖 Citation

> Your Name. SPLINTR: A tool to predict coding disruption and NMD from alternative splicing. GitHub. 2024.
