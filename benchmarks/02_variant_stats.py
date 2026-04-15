"""
Benchmark 2 — Variant quality statistics (works on ANY 10 samples).

For each sample's PASS VCF, uses bcftools stats to calculate:
  - Total variants, SNPs, INDELs
  - Ts/Tv ratio  (expected 2.0–3.3 for exomes)
  - Het / Hom ratio  (expected ~2.0 for germline)
  - dbSNP concordance (% of variants in dbSNP)
  - Mean QUAL score
  - Variant counts per chromosome

Generates:
  - benchmarks/results/variant_stats.csv
  - benchmarks/results/variant_counts.png    (bar chart SNP/INDEL per sample)
    - benchmarks/results/tstv_hethom.png       (Ts/Tv and Het/Hom QC plot)

Usage
-----
python benchmarks/02_variant_stats.py \
    --vcf-dir   results/VCF \
    --dbsnp     /path/to/dbsnp.vcf.gz

Requirements
------------
bcftools must be on PATH
pip install pandas matplotlib seaborn
"""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import seaborn as sns


def run_bcftools_stats(vcf: Path, dbsnp: Path | None) -> dict:
    """Run bcftools stats and parse key metrics."""
    cmd = ["bcftools", "stats"]
    if dbsnp:
        cmd += ["--dbsnp", str(dbsnp)]
    cmd.append(str(vcf))

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"bcftools stats failed for {vcf.name}:\n{result.stderr}")

    out = result.stdout
    metrics: dict = {}

    def val(pattern: str) -> str:
        m = re.search(pattern, out, re.MULTILINE)
        return m.group(1).strip() if m else "0"

    # Summary numbers block (SN lines)
    metrics["total_variants"] = int(val(r"^SN\t0\tnumber of records:\t(\d+)"))
    metrics["snps"]           = int(val(r"^SN\t0\tnumber of SNPs:\t(\d+)"))
    metrics["indels"]         = int(val(r"^SN\t0\tnumber of indels:\t(\d+)"))
    metrics["mnps"]           = int(val(r"^SN\t0\tnumber of MNPs:\t(\d+)"))
    metrics["multiallelic"]   = int(val(r"^SN\t0\tnumber of multiallelic sites:\t(\d+)"))

    # Ts/Tv
    tstv_m = re.search(r"^TSTV\t0\t[\d.]+\t[\d.]+\t([\d.]+)", out, re.MULTILINE)
    metrics["tstv"] = float(tstv_m.group(1)) if tstv_m else 0.0

    # Het/Hom
    het_m  = re.search(r"^SN\t0\tnumber of het SNPs:\t(\d+)", out, re.MULTILINE)
    hom_m  = re.search(r"^SN\t0\tnumber of hom SNPs:\t(\d+)", out, re.MULTILINE)
    n_het = int(het_m.group(1)) if het_m else 0
    n_hom = int(hom_m.group(1)) if hom_m else 1
    metrics["het_snps"]  = n_het
    metrics["hom_snps"]  = n_hom
    metrics["het_hom_ratio"] = round(n_het / n_hom, 3) if n_hom else 0.0

    # dbSNP concordance (only available when --dbsnp is passed)
    dbsnp_m = re.search(r"^SN\t0\tnumber of SNPs in dbSNP:\t(\d+)", out, re.MULTILINE)
    if dbsnp_m and metrics["snps"] > 0:
        metrics["dbsnp_concordance_pct"] = round(
            int(dbsnp_m.group(1)) / metrics["snps"] * 100, 2
        )
    else:
        metrics["dbsnp_concordance_pct"] = None

    return metrics


def collect_all(vcf_dir: Path, dbsnp: Path | None) -> pd.DataFrame:
    records = []
    vcfs = sorted(vcf_dir.glob("*_PASS.vcf"))
    if not vcfs:
        raise SystemExit(
            f"No *_PASS.vcf files found in {vcf_dir}.\n"
            "Run the pipeline first: exomeflow run ..."
        )
    for vcf in vcfs:
        sample = vcf.stem.replace("_PASS", "")
        print(f"  Processing {sample} ...", end="", flush=True)
        try:
            m = run_bcftools_stats(vcf, dbsnp)
            m["sample"] = sample
            records.append(m)
            print(f" SNPs={m['snps']:,}  INDELs={m['indels']:,}  Ts/Tv={m['tstv']:.2f}")
        except RuntimeError as e:
            print(f" ERROR: {e}")
    return pd.DataFrame(records).set_index("sample")


def plot_counts(df: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(14, 5))
    x = range(len(df))
    width = 0.35
    ax.bar([i - width / 2 for i in x], df["snps"],   width, label="SNPs",   color="#2196F3")
    ax.bar([i + width / 2 for i in x], df["indels"], width, label="INDELs", color="#FF9800")
    ax.set_xticks(list(x))
    ax.set_xticklabels(df.index, rotation=30, ha="right")
    ax.set_ylabel("Variant Count")
    ax.set_title("ExomeFlow — PASS Variant Counts per Sample")
    ax.legend()
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))
    plt.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"  Saved: {out}")


def plot_qc(df: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Ts/Tv ratio — highlight acceptable range
    axes[0].bar(df.index, df["tstv"], color="#4CAF50", edgecolor="white")
    axes[0].axhspan(2.0, 3.3, alpha=0.15, color="green", label="Expected range (2.0–3.3)")
    axes[0].set_title("Ts/Tv Ratio")
    axes[0].set_ylabel("Ts/Tv")
    axes[0].set_xticklabels(df.index, rotation=30, ha="right")
    axes[0].legend(fontsize=8)
    for bar, val in zip(axes[0].patches, df["tstv"]):
        axes[0].text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                     f"{val:.2f}", ha="center", va="bottom", fontsize=7)

    # Het/Hom ratio — highlight acceptable range
    axes[1].bar(df.index, df["het_hom_ratio"], color="#9C27B0", edgecolor="white")
    axes[1].axhspan(1.5, 2.5, alpha=0.15, color="purple", label="Expected range (1.5–2.5)")
    axes[1].set_title("Het/Hom SNP Ratio")
    axes[1].set_ylabel("Het/Hom Ratio")
    axes[1].set_xticklabels(df.index, rotation=30, ha="right")
    axes[1].legend(fontsize=8)
    for bar, val in zip(axes[1].patches, df["het_hom_ratio"]):
        axes[1].text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                     f"{val:.2f}", ha="center", va="bottom", fontsize=7)

    fig.suptitle("ExomeFlow — Variant Quality Control Metrics", fontsize=13)
    plt.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"  Saved: {out}")


def flag_outliers(df: pd.DataFrame) -> None:
    """Print QC warnings for samples outside expected ranges."""
    print("\n  QC Flags:")
    flagged = False
    for sample, row in df.iterrows():
        flags = []
        if row["tstv"] < 2.0 or row["tstv"] > 3.3:
            flags.append(f"Ts/Tv={row['tstv']:.2f} (expected 2.0–3.3)")
        if row["het_hom_ratio"] < 1.5 or row["het_hom_ratio"] > 2.5:
            flags.append(f"Het/Hom={row['het_hom_ratio']:.2f} (expected 1.5–2.5)")
        if row["snps"] < 20000:
            flags.append(f"Low SNP count={row['snps']:,}")
        if flags:
            print(f"    WARN  {sample}: {' | '.join(flags)}")
            flagged = True
    if not flagged:
        print("    All samples within expected ranges.")


def main() -> None:
    ap = argparse.ArgumentParser(description="ExomeFlow variant quality stats")
    ap.add_argument("--vcf-dir",  required=True, type=Path,
                    help="Directory containing *_PASS.vcf files (e.g. results/VCF)")
    ap.add_argument("--dbsnp",    default=None,  type=Path,
                    help="dbSNP VCF for concordance calculation (optional but recommended)")
    ap.add_argument("--out-dir",  default="benchmarks/results", type=Path)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nRunning bcftools stats on VCFs in: {args.vcf_dir}\n")
    df = collect_all(args.vcf_dir, args.dbsnp)

    # Print summary table
    cols = ["total_variants", "snps", "indels", "tstv", "het_hom_ratio", "dbsnp_concordance_pct"]
    cols = [c for c in cols if c in df.columns]
    print("\n" + df[cols].to_string())
    flag_outliers(df)

    # Save CSV
    csv_path = args.out_dir / "variant_stats.csv"
    df.to_csv(csv_path)
    print(f"\n  Saved: {csv_path}")

    # Figures
    plot_counts(df, args.out_dir / "variant_counts.png")
    plot_qc(df,     args.out_dir / "tstv_hethom.png")
    print("\nDone.")


if __name__ == "__main__":
    main()
