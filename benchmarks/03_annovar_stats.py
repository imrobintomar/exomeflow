"""
Benchmark 3 — Functional annotation statistics from ANNOVAR output.

For each sample's hg38_multianno.txt file, calculates:
  - Variant functional consequence distribution (exonic / intronic / splicing …)
  - Exonic function breakdown (synonymous / nonsynonymous / stopgain …)
  - ClinVar classification distribution
  - gnomAD allele frequency distribution (novel vs rare vs common)
  - Novel variants (not in dbSNP / avSNP150)
  - Cross-sample comparison of pathogenic/likely-pathogenic variant counts

Generates:
  - benchmarks/results/annovar_stats.csv
  - benchmarks/results/func_distribution.png    (stacked bar per sample)
  - benchmarks/results/clinvar_distribution.png (ClinVar classifications)
  - benchmarks/results/gnomad_af.png            (AF spectrum per sample)
  - benchmarks/results/novel_vs_known.png       (novel/known/rare breakdown)

Usage
-----
python benchmarks/03_annovar_stats.py \
    --vcf-dir results/VCF \
    --out-dir benchmarks/results

Requirements
------------
pip install pandas matplotlib seaborn
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns


# ANNOVAR column names (hg38_multianno.txt header)
COL_FUNC        = "Func.refGene"
COL_EXFUNC      = "ExonicFunc.refGene"
COL_GENE        = "Gene.refGene"
COL_CLINVAR     = "CLNSIG"
COL_GNOMAD_EX   = "gnomAD_exome_ALL"
COL_GNOMAD_GEN  = "gnomAD_genome_ALL"
COL_AVSNP       = "avsnp150"
COL_CADD        = "CADD_phred"

# Functional consequence colour map
FUNC_COLORS = {
    "exonic":           "#E53935",
    "splicing":         "#FB8C00",
    "exonic;splicing":  "#FF7043",
    "ncRNA_exonic":     "#8E24AA",
    "ncRNA_splicing":   "#AB47BC",
    "UTR5":             "#1E88E5",
    "UTR3":             "#42A5F5",
    "intronic":         "#78909C",
    "upstream":         "#A5D6A7",
    "downstream":       "#C8E6C9",
    "intergenic":       "#EEEEEE",
}

CLINVAR_COLORS = {
    "Pathogenic":                "#C62828",
    "Likely_pathogenic":         "#EF5350",
    "Pathogenic/Likely_pathogenic": "#E53935",
    "Uncertain_significance":    "#FB8C00",
    "Likely_benign":             "#66BB6A",
    "Benign":                    "#2E7D32",
    "Benign/Likely_benign":      "#388E3C",
    "not_provided":              "#BDBDBD",
    "other":                     "#E0E0E0",
}

EXFUNC_COLORS = {
    "nonsynonymous SNV":  "#C62828",
    "stopgain":           "#B71C1C",
    "stoploss":           "#E53935",
    "startloss":          "#FF5252",
    "frameshift insertion":  "#FF6D00",
    "frameshift deletion":   "#FF9100",
    "nonframeshift insertion": "#FFD600",
    "nonframeshift deletion":  "#FFFF00",
    "synonymous SNV":     "#81C784",
    "unknown":            "#BDBDBD",
}


def load_multianno(path: Path) -> pd.DataFrame:
    """Load ANNOVAR multianno.txt, normalise column names."""
    df = pd.read_csv(path, sep="\t", low_memory=False, dtype=str)
    df.columns = df.columns.str.strip()
    # Replace '.' with NaN
    df.replace(".", np.nan, inplace=True)
    return df


def func_counts(df: pd.DataFrame) -> pd.Series:
    if COL_FUNC not in df.columns:
        return pd.Series(dtype=int)
    counts = df[COL_FUNC].fillna("unknown").str.split(";").explode().value_counts()
    return counts


def exfunc_counts(df: pd.DataFrame) -> pd.Series:
    if COL_EXFUNC not in df.columns:
        return pd.Series(dtype=int)
    exonic = df[df[COL_FUNC].fillna("").str.contains("exonic", case=False)]
    return exonic[COL_EXFUNC].fillna("unknown").value_counts()


def clinvar_counts(df: pd.DataFrame) -> pd.Series:
    for col in [COL_CLINVAR, "CLNALLELEID", "CLNREVSTAT"]:
        if col == COL_CLINVAR and col in df.columns:
            raw = df[col].dropna()
            raw = raw[raw != ""]
            # Normalise compound labels like "Pathogenic,Likely_pathogenic"
            series = (
                raw.str.replace(r"[,|]", "/", regex=True)
                   .str.strip()
                   .value_counts()
            )
            # Map uncommon labels to "other"
            known = set(CLINVAR_COLORS.keys())
            series.index = [i if i in known else "other" for i in series.index]
            return series.groupby(series.index).sum().sort_values(ascending=False)
    return pd.Series(dtype=int)


def gnomad_af_bins(df: pd.DataFrame) -> pd.Series:
    """Bin variants by gnomAD AF: novel / ultra-rare / rare / low-freq / common."""
    for col in [COL_GNOMAD_EX, COL_GNOMAD_GEN]:
        if col in df.columns:
            af = pd.to_numeric(df[col], errors="coerce")
            bins   = [-1, 0, 0.0001, 0.001, 0.01, 0.05, 1.01]
            labels = ["Novel", "Ultra-rare\n(<0.01%)", "Rare\n(0.01–0.1%)",
                      "Low-freq\n(0.1–1%)", "Common\n(1–5%)", "Very common\n(>5%)"]
            return pd.cut(af.fillna(-1), bins=bins, labels=labels).value_counts().sort_index()
    return pd.Series(dtype=int)


def novel_known(df: pd.DataFrame) -> dict:
    """Count novel (no dbSNP ID) vs known variants."""
    if COL_AVSNP not in df.columns:
        return {}
    total  = len(df)
    known  = df[COL_AVSNP].notna().sum()
    novel  = total - known
    return {"Novel": int(novel), "Known (dbSNP)": int(known)}


def collect_all(vcf_dir: Path) -> dict[str, pd.DataFrame]:
    """Load all multianno.txt files. Returns {sample: dataframe}."""
    files = sorted(vcf_dir.glob("*.annovar.hg38_multianno.txt"))
    if not files:
        raise SystemExit(
            f"No *.annovar.hg38_multianno.txt files found in {vcf_dir}.\n"
            "Run the pipeline first: exomeflow run ..."
        )
    data = {}
    for f in files:
        sample = re.sub(r"\.annovar\.hg38_multianno\.txt$", "", f.name)
        print(f"  Loading {sample} ({f.stat().st_size // 1024:,} KB) ...", flush=True)
        data[sample] = load_multianno(f)
    return data


def plot_func_distribution(sample_data: dict[str, pd.DataFrame], out: Path) -> None:
    rows = {}
    for sample, df in sample_data.items():
        rows[sample] = func_counts(df)
    combined = pd.DataFrame(rows).fillna(0).astype(int).T

    # Keep top-N categories
    top_cats = combined.sum().nlargest(10).index
    combined = combined[top_cats]

    colors = [FUNC_COLORS.get(c, "#BDBDBD") for c in combined.columns]
    fig, ax = plt.subplots(figsize=(14, 6))
    combined.plot.bar(stacked=True, ax=ax, color=colors, edgecolor="white", width=0.75)
    ax.set_title("Variant Functional Consequence Distribution per Sample")
    ax.set_ylabel("Variant Count")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(loc="upper right", fontsize=7, ncol=2)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))
    plt.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"  Saved: {out}")


def plot_clinvar(sample_data: dict[str, pd.DataFrame], out: Path) -> None:
    rows = {}
    for sample, df in sample_data.items():
        c = clinvar_counts(df)
        if not c.empty:
            rows[sample] = c
    if not rows:
        print("  No ClinVar data found — skipping clinvar_distribution.png")
        return

    combined = pd.DataFrame(rows).fillna(0).astype(int).T
    path_cols = [c for c in ["Pathogenic", "Likely_pathogenic",
                              "Pathogenic/Likely_pathogenic"] if c in combined.columns]
    benign_cols = [c for c in ["Benign", "Likely_benign",
                                "Benign/Likely_benign"] if c in combined.columns]
    other_cols  = [c for c in combined.columns if c not in path_cols + benign_cols]
    order = path_cols + benign_cols + other_cols
    combined = combined[[c for c in order if c in combined.columns]]

    colors = [CLINVAR_COLORS.get(c, "#E0E0E0") for c in combined.columns]
    fig, ax = plt.subplots(figsize=(14, 6))
    combined.plot.bar(stacked=True, ax=ax, color=colors, edgecolor="white", width=0.75)
    ax.set_title("ClinVar Classification Distribution per Sample")
    ax.set_ylabel("Variant Count")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(loc="upper right", fontsize=7, ncol=2)
    plt.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"  Saved: {out}")


def plot_gnomad_af(sample_data: dict[str, pd.DataFrame], out: Path) -> None:
    rows = {}
    for sample, df in sample_data.items():
        rows[sample] = gnomad_af_bins(df)
    combined = pd.DataFrame(rows).fillna(0).astype(int).T

    if combined.empty:
        print("  No gnomAD AF data found — skipping gnomad_af.png")
        return

    palette = ["#C62828", "#EF5350", "#FB8C00", "#FFF176", "#AED581", "#388E3C"]
    fig, ax = plt.subplots(figsize=(14, 6))
    combined.plot.bar(stacked=True, ax=ax,
                      color=palette[:len(combined.columns)],
                      edgecolor="white", width=0.75)
    ax.set_title("gnomAD Allele Frequency Spectrum per Sample")
    ax.set_ylabel("Variant Count")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(loc="upper right", fontsize=7)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))
    plt.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"  Saved: {out}")


def plot_novel_known(sample_data: dict[str, pd.DataFrame], out: Path) -> None:
    rows = {}
    for sample, df in sample_data.items():
        nk = novel_known(df)
        if nk:
            rows[sample] = nk
    if not rows:
        print("  No avSNP150 data — skipping novel_vs_known.png")
        return

    combined = pd.DataFrame(rows).T.fillna(0).astype(int)
    fig, ax = plt.subplots(figsize=(12, 5))
    combined.plot.bar(stacked=True, ax=ax,
                      color=["#E53935", "#1E88E5"],
                      edgecolor="white", width=0.7)
    ax.set_title("Novel vs Known (dbSNP) Variants per Sample")
    ax.set_ylabel("Variant Count")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(fontsize=9)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))

    # Add % novel labels
    for i, (sample, row) in enumerate(combined.iterrows()):
        total = row.sum()
        pct   = row.get("Novel", 0) / total * 100 if total > 0 else 0
        ax.text(i, total + total * 0.01, f"{pct:.1f}%\nnovel",
                ha="center", va="bottom", fontsize=7)

    plt.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"  Saved: {out}")


def build_summary_csv(sample_data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    records = []
    for sample, df in sample_data.items():
        row = {"sample": sample, "total_variants": len(df)}

        # Functional counts
        fc = func_counts(df)
        row["exonic"]    = int(fc.get("exonic", 0))
        row["splicing"]  = int(fc.get("splicing", 0))
        row["intronic"]  = int(fc.get("intronic", 0))

        # Exonic function
        ef = exfunc_counts(df)
        row["nonsynonymous"] = int(ef.get("nonsynonymous SNV", 0))
        row["synonymous"]    = int(ef.get("synonymous SNV", 0))
        row["stopgain"]      = int(ef.get("stopgain", 0))
        row["frameshift"]    = int(
            ef.filter(like="frameshift").sum() if not ef.empty else 0
        )

        # ClinVar pathogenic
        cv = clinvar_counts(df)
        row["pathogenic"] = int(
            sum(cv.get(k, 0) for k in
                ["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"])
        )

        # gnomAD novel
        af_bins = gnomad_af_bins(df)
        row["novel_variants"] = int(af_bins.get("Novel", 0))

        # dbSNP novel
        nk = novel_known(df)
        row["novel_dbsnp"] = nk.get("Novel", "N/A")

        records.append(row)
    return pd.DataFrame(records).set_index("sample")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="ExomeFlow ANNOVAR functional annotation statistics"
    )
    ap.add_argument("--vcf-dir",  required=True, type=Path,
                    help="Directory containing *.annovar.hg38_multianno.txt files")
    ap.add_argument("--out-dir",  default="benchmarks/results", type=Path)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nLoading ANNOVAR results from: {args.vcf_dir}\n")
    sample_data = collect_all(args.vcf_dir)

    print(f"\nBuilding summary table ...")
    summary = build_summary_csv(sample_data)
    print("\n" + summary.to_string())

    csv_path = args.out_dir / "annovar_stats.csv"
    summary.to_csv(csv_path)
    print(f"\n  Saved: {csv_path}")

    print("\nGenerating figures ...")
    plot_func_distribution(sample_data, args.out_dir / "func_distribution.png")
    plot_clinvar(sample_data,           args.out_dir / "clinvar_distribution.png")
    plot_gnomad_af(sample_data,         args.out_dir / "gnomad_af.png")
    plot_novel_known(sample_data,       args.out_dir / "novel_vs_known.png")

    print("\nDone.")


if __name__ == "__main__":
    main()
