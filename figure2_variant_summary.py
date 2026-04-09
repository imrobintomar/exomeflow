"""
Figure 2 — Variant Summary Bar Chart for ExomeFlow paper.

Reads directly from VCF output files produced by the pipeline and
generates a publication-quality grouped bar chart.

Usage:
    python figure2_variant_summary.py

Output:
    figure2_variant_summary.png  (300 dpi, for publication)
    figure2_variant_summary.svg  (vector, for editing in Inkscape)
    variant_counts.csv           (data table for supplementary)
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

# ── CONFIG ────────────────────────────────────────────────────────────────────
VCF_DIR = Path("/media/drprabudh/m2/ARM-PRJNA1344179/fastq/VCF")
OUTPUT_DIR = Path("/media/drprabudh/m2/Exomeflow")

# Short labels for the x-axis (strip the SRR prefix for readability)
LABEL_PREFIX = "SRR3679"   # everything after this becomes the tick label


# ── VARIANT COUNTING ──────────────────────────────────────────────────────────

def count_vcf_lines(vcf: Path) -> int:
    """Count non-header variant lines in a VCF file."""
    if not vcf.exists():
        return 0
    count = 0
    with open(vcf, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("#"):
                count += 1
    return count


def classify_variants(vcf: Path) -> tuple[int, int]:
    """
    Return (snp_count, indel_count) by inspecting REF/ALT columns.

    A variant is classified as SNP if both REF and ALT are single
    nucleotides (after splitting multi-allelic ALT on comma).
    Everything else is counted as INDEL.
    """
    if not vcf.exists():
        return 0, 0

    snps = 0
    indels = 0

    with open(vcf, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            ref = parts[3].strip()
            alts = parts[4].strip().split(",")
            is_snp = len(ref) == 1 and all(len(a) == 1 and a != "*" for a in alts)
            if is_snp:
                snps += 1
            else:
                indels += 1

    return snps, indels


# ── BUILD DATA TABLE ──────────────────────────────────────────────────────────

def build_table(vcf_dir: Path) -> pd.DataFrame:
    raw_vcfs = sorted(vcf_dir.glob("SRR*.vcf"))
    # Exclude PASS / annovar files
    raw_vcfs = [v for v in raw_vcfs if "_PASS" not in v.name and "annovar" not in v.name]

    records = []
    for raw in raw_vcfs:
        sample = raw.stem                              # e.g. SRR36790057
        pass_vcf = vcf_dir / f"{sample}_PASS.vcf"

        total         = count_vcf_lines(raw)
        snps, indels  = classify_variants(raw)
        passed        = count_vcf_lines(pass_vcf)
        filtered      = total - passed

        records.append({
            "sample":   sample,
            "label":    sample.replace(LABEL_PREFIX, "005"),  # short x-axis label
            "total":    total,
            "snps":     snps,
            "indels":   indels,
            "pass":     passed,
            "filtered": filtered,
        })

    return pd.DataFrame(records)


# ── PLOT ──────────────────────────────────────────────────────────────────────

def plot(df: pd.DataFrame, output_dir: Path) -> None:
    # ── style
    plt.rcParams.update({
        "font.family":     "Arial",
        "font.size":       11,
        "axes.linewidth":  1.2,
        "axes.spines.top":    False,
        "axes.spines.right":  False,
        "xtick.direction": "out",
        "ytick.direction": "out",
    })

    n       = len(df)
    x       = range(n)
    width   = 0.20          # bar width
    offsets = [-1.5, -0.5, 0.5, 1.5]   # 4 groups per sample

    # Colour palette (colour-blind friendly)
    colours = {
        "total":    "#4878CF",   # blue
        "snps":     "#6ACC65",   # green
        "indels":   "#D65F5F",   # red
        "pass":     "#B47CC7",   # purple
    }

    fig, ax = plt.subplots(figsize=(12, 6))

    for i, (key, offset, label_text) in enumerate([
        ("total",   offsets[0], "Total Variants"),
        ("snps",    offsets[1], "SNPs"),
        ("indels",  offsets[2], "INDELs"),
        ("pass",    offsets[3], "PASS Variants"),
    ]):
        positions = [xi + offset * width for xi in x]
        bars = ax.bar(
            positions,
            df[key],
            width=width,
            color=colours[key],
            edgecolor="white",
            linewidth=0.6,
            label=label_text,
            zorder=3,
        )

        # Value labels on top of each bar
        for bar in bars:
            h = bar.get_height()
            if h > 0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    h + df[key].max() * 0.008,
                    f"{int(h):,}",
                    ha="center", va="bottom",
                    fontsize=7.5, rotation=90,
                    color="#333333",
                )

    # ── axes
    ax.set_xticks(list(x))
    ax.set_xticklabels(df["sample"].tolist(), rotation=30, ha="right", fontsize=10)
    ax.set_ylabel("Number of Variants", fontsize=12, labelpad=10)
    ax.set_xlabel("Sample (SRA Accession)", fontsize=12, labelpad=10)
    ax.set_title(
        "Variant Summary per Sample after ExomeFlow Analysis",
        fontsize=13, fontweight="bold", pad=14,
    )

    ax.yaxis.grid(True, linestyle="--", alpha=0.5, zorder=0)
    ax.set_axisbelow(True)

    # ── legend
    ax.legend(
        frameon=False,
        fontsize=10,
        loc="upper right",
        ncol=2,
    )

    # ── filtered annotation (small text below x-axis)
    for xi, row in zip(x, df.itertuples()):
        ax.text(
            xi, -df["total"].max() * 0.07,
            f"Filtered:\n{row.filtered:,}",
            ha="center", va="top",
            fontsize=7, color="#888888",
        )

    plt.tight_layout()

    png_path = output_dir / "figure2_variant_summary.png"
    svg_path = output_dir / "figure2_variant_summary.svg"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    print(f"Saved: {png_path}")
    print(f"Saved: {svg_path}")
    plt.show()


# ── MAIN ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("Counting variants from VCF files ...")
    df = build_table(VCF_DIR)

    # Save data table
    csv_path = OUTPUT_DIR / "variant_counts.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nVariant counts:\n{df[['sample','total','snps','indels','pass','filtered']].to_string(index=False)}")
    print(f"\nSaved data table: {csv_path}")

    print("\nGenerating Figure 2 ...")
    plot(df, OUTPUT_DIR)
