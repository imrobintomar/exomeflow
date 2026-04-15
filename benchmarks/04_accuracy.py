"""
Benchmark 3 — Accuracy benchmarking using hap.py (requires NA12878 / GIAB sample).

For each sample VCF, compares against the GIAB truth set using hap.py and
calculates sensitivity, precision, and F1 score for SNPs and INDELs separately.

Generates:
  - benchmarks/results/accuracy_table.csv
  - benchmarks/results/accuracy_snp_indel.png   (precision/recall/F1 bar chart)
  - benchmarks/results/roc_curve.png             (ROC-style plot)

Prerequisites
-------------
1. hap.py installed:
       pip install hap.py
   or:  conda install -c bioconda hap.py

2. GIAB truth VCF + high-confidence BED for hg38:
       wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
       wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed

3. Your reference genome FASTA (hg38.fa)

Usage
-----
# Single sample (NA12878):
python benchmarks/03_accuracy.py \
    --query-vcf   results/VCF/NA12878_PASS.vcf \
    --truth-vcf   /data/refs/HG001_GRCh38_v4.2.1_benchmark.vcf.gz \
    --truth-bed   /data/refs/HG001_GRCh38_v4.2.1_benchmark.bed \
    --reference   /data/refs/hg38.fa \
    --sample-name NA12878

# Multiple samples (all matching a pattern):
python benchmarks/03_accuracy.py \
    --vcf-dir     results/VCF \
    --truth-vcf   /data/refs/HG001_GRCh38_v4.2.1_benchmark.vcf.gz \
    --truth-bed   /data/refs/HG001_GRCh38_v4.2.1_benchmark.bed \
    --reference   /data/refs/hg38.fa
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd


GIAB_TRUTH_VCF = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
)
GIAB_TRUTH_BED = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
)


def check_hap_py() -> str:
    """Return path to hap.py executable, or raise."""
    for candidate in ["hap.py", "happy"]:
        if shutil.which(candidate):
            return candidate
    raise SystemExit(
        "hap.py not found on PATH.\n"
        "Install: conda install -c bioconda hap.py\n"
        "Or:      pip install hap.py"
    )


def run_hap_py(
    query_vcf: Path,
    truth_vcf: Path,
    truth_bed: Path,
    reference: Path,
    out_prefix: Path,
    hap_py: str,
) -> dict:
    """Run hap.py and parse the extended CSV output."""
    cmd = [
        hap_py,
        str(truth_vcf),
        str(query_vcf),
        "-f", str(truth_bed),
        "-r", str(reference),
        "-o", str(out_prefix),
        "--engine", "vcfeval",
        "--pass-only",
        "--no-roc",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"hap.py failed:\n{result.stderr[-2000:]}")

    # Parse extended CSV (hap.py writes <prefix>.extended.csv)
    csv_path = Path(str(out_prefix) + ".extended.csv")
    if not csv_path.exists():
        # Fall back to summary.csv
        csv_path = Path(str(out_prefix) + ".summary.csv")
    df = pd.read_csv(csv_path)

    metrics = {}
    for _, row in df.iterrows():
        vtype = row.get("Type", "")
        filt  = row.get("Filter", "PASS")
        if filt != "PASS":
            continue
        if vtype in ("SNP", "INDEL"):
            metrics[vtype] = {
                "precision": round(float(row.get("METRIC.Precision", 0)), 4),
                "recall":    round(float(row.get("METRIC.Recall",    0)), 4),
                "f1":        round(float(row.get("METRIC.F1_Score",  0)), 4),
                "tp":        int(row.get("TRUTH.TP", 0)),
                "fn":        int(row.get("TRUTH.FN", 0)),
                "fp":        int(row.get("QUERY.FP", 0)),
            }
    return metrics


def collect_all(
    vcf_dir: Path,
    truth_vcf: Path,
    truth_bed: Path,
    reference: Path,
    out_dir: Path,
    hap_py: str,
) -> pd.DataFrame:
    records = []
    for vcf in sorted(vcf_dir.glob("*_PASS.vcf")):
        sample = vcf.stem.replace("_PASS", "")
        print(f"  Running hap.py on {sample} ...", flush=True)
        prefix = out_dir / "happy" / sample
        prefix.parent.mkdir(parents=True, exist_ok=True)
        try:
            m = run_hap_py(vcf, truth_vcf, truth_bed, reference, prefix, hap_py)
            for vtype, vals in m.items():
                row = {"sample": sample, "type": vtype}
                row.update(vals)
                records.append(row)
            snp = m.get("SNP", {})
            ind = m.get("INDEL", {})
            print(
                f"    SNP  : Precision={snp.get('precision', 0):.3f}  "
                f"Recall={snp.get('recall', 0):.3f}  F1={snp.get('f1', 0):.3f}"
            )
            print(
                f"    INDEL: Precision={ind.get('precision', 0):.3f}  "
                f"Recall={ind.get('recall', 0):.3f}  F1={ind.get('f1', 0):.3f}"
            )
        except RuntimeError as e:
            print(f"    ERROR: {e}")
    if not records:
        raise SystemExit("No results — check hap.py errors above.")
    return pd.DataFrame(records)


def plot_accuracy(df: pd.DataFrame, out: Path) -> None:
    metrics = ["precision", "recall", "f1"]
    vtypes  = ["SNP", "INDEL"]
    colors  = {"precision": "#2196F3", "recall": "#4CAF50", "f1": "#FF5722"}

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=False)
    for ax, vtype in zip(axes, vtypes):
        sub = df[df["type"] == vtype].set_index("sample")
        x   = range(len(sub))
        width = 0.25
        for i, metric in enumerate(metrics):
            if metric in sub.columns:
                ax.bar(
                    [xi + i * width for xi in x],
                    sub[metric],
                    width, label=metric.capitalize(),
                    color=colors[metric], edgecolor="white",
                )
        ax.set_xticks([xi + width for xi in x])
        ax.set_xticklabels(sub.index, rotation=30, ha="right", fontsize=8)
        ax.set_ylim(0, 1.05)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0])
        ax.set_ylabel("Score")
        ax.set_title(f"{vtype} Accuracy")
        ax.legend(fontsize=8)
        ax.axhline(0.95, ls="--", lw=0.8, color="gray", label="0.95 target")

    fig.suptitle("ExomeFlow — Variant Calling Accuracy (vs GIAB NA12878)", fontsize=13)
    plt.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"  Saved: {out}")


def plot_precision_recall(df: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 7))
    colors  = {"SNP": "#2196F3", "INDEL": "#FF9800"}
    markers = {"SNP": "o", "INDEL": "s"}

    for vtype in ["SNP", "INDEL"]:
        sub = df[df["type"] == vtype]
        ax.scatter(sub["recall"], sub["precision"],
                   c=colors[vtype], marker=markers[vtype],
                   s=80, label=vtype, zorder=3)
        for _, row in sub.iterrows():
            ax.annotate(row["sample"],
                        (row["recall"], row["precision"]),
                        textcoords="offset points", xytext=(5, 3),
                        fontsize=6, color="gray")

    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Recall (Sensitivity)")
    ax.set_ylabel("Precision (PPV)")
    ax.set_title("Precision vs Recall — ExomeFlow vs GIAB Truth")
    ax.axhline(0.95, ls="--", lw=0.7, color="gray")
    ax.axvline(0.95, ls="--", lw=0.7, color="gray")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"  Saved: {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description="ExomeFlow accuracy benchmark (hap.py)")
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--query-vcf", type=Path,
                     help="Single query VCF (PASS variants)")
    grp.add_argument("--vcf-dir",   type=Path,
                     help="Directory containing *_PASS.vcf files")
    ap.add_argument("--truth-vcf",   required=True, type=Path,
                    help="GIAB truth VCF (bgzipped .vcf.gz)")
    ap.add_argument("--truth-bed",   required=True, type=Path,
                    help="GIAB high-confidence regions BED")
    ap.add_argument("--reference",   required=True, type=Path,
                    help="Reference genome FASTA (hg38.fa)")
    ap.add_argument("--sample-name", default=None,
                    help="Sample name (used with --query-vcf)")
    ap.add_argument("--out-dir",     default="benchmarks/results", type=Path)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    hap_py = check_hap_py()
    print(f"\nUsing hap.py: {hap_py}")
    print(f"Truth VCF   : {args.truth_vcf}")
    print(f"Truth BED   : {args.truth_bed}\n")

    if args.query_vcf:
        name = args.sample_name or args.query_vcf.stem.replace("_PASS", "")
        prefix = args.out_dir / "happy" / name
        prefix.parent.mkdir(parents=True, exist_ok=True)
        m = run_hap_py(args.query_vcf, args.truth_vcf,
                       args.truth_bed, args.reference, prefix, hap_py)
        records = []
        for vtype, vals in m.items():
            row = {"sample": name, "type": vtype}
            row.update(vals)
            records.append(row)
        df = pd.DataFrame(records)
    else:
        df = collect_all(args.vcf_dir, args.truth_vcf,
                         args.truth_bed, args.reference, args.out_dir, hap_py)

    # Print summary
    print("\n" + df.to_string(index=False))

    # Save CSV
    csv_path = args.out_dir / "accuracy_table.csv"
    df.to_csv(csv_path, index=False, float_format="%.4f")
    print(f"\n  Saved: {csv_path}")

    # Figures
    plot_accuracy(df, args.out_dir / "accuracy_snp_indel.png")
    plot_precision_recall(df, args.out_dir / "precision_recall.png")
    print("\nDone.")

    # Final summary line
    for vtype in ["SNP", "INDEL"]:
        sub = df[df["type"] == vtype]
        if not sub.empty:
            print(
                f"  Mean {vtype:5s}  Precision={sub['precision'].mean():.3f}  "
                f"Recall={sub['recall'].mean():.3f}  F1={sub['f1'].mean():.3f}"
            )


if __name__ == "__main__":
    main()
