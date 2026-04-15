"""
Benchmark 1 — Performance (runtime per step, per sample).

Parses ExomeFlow log files to extract wall-clock time for every pipeline
step, then generates:
  - benchmarks/results/runtime_table.csv
  - benchmarks/results/runtime_heatmap.png
  - benchmarks/results/runtime_bars.png

Usage
-----
python benchmarks/01_performance.py --log-dir results/logs

Requirements
------------
pip install pandas matplotlib seaborn
"""

from __future__ import annotations

import argparse
import re
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import seaborn as sns

# Step names as they appear in log messages
STEP_PATTERNS = {
    "fastp":          r"fastp.*completed",
    "BWA MEM":        r"BWA MEM.*completed",
    "SortSam":        r"SortSam.*completed|BAM sort.*completed",
    "flagstat":       r"flagstat.*completed",
    "MarkDuplicates": r"MarkDuplicates.*completed|duplicate.*completed",
    "BuildBamIndex":  r"BuildBamIndex.*completed|BAM index.*completed",
    "BQSR":           r"BQSR.*completed",
    "HaplotypeCaller":r"HaplotypeCaller.*completed|variant calling.*completed",
    "Filtering":      r"[Ff]iltering.*completed",
    "ANNOVAR":        r"ANNOVAR.*completed|annotation.*completed",
}

START_PATTERNS = {
    "fastp":           r"Running fastp",
    "BWA MEM":         r"Running BWA MEM|Aligning",
    "SortSam":         r"SortSam|Sorting BAM",
    "flagstat":        r"flagstat",
    "MarkDuplicates":  r"MarkDuplicates|Marking duplicates",
    "BuildBamIndex":   r"BuildBamIndex|Indexing BAM",
    "BQSR":            r"BaseRecalibrator|Running BQSR",
    "HaplotypeCaller": r"HaplotypeCaller|Running.*variant call",
    "Filtering":       r"Separating SNPs|Running.*filter",
    "ANNOVAR":         r"Running ANNOVAR",
}

TIMESTAMP_FMT = "%Y-%m-%d %H:%M:%S"
LOG_RE = re.compile(r"\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\].*\[(\w+)\]\s+(.*)")


def parse_log(log_path: Path) -> dict[str, float]:
    """Return {step: duration_seconds} for one sample log file."""
    entries: list[tuple[datetime, str]] = []
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = LOG_RE.match(line)
        if m:
            ts = datetime.strptime(m.group(1), TIMESTAMP_FMT)
            msg = m.group(3)
            entries.append((ts, msg))

    durations: dict[str, float] = {}
    for step, start_pat in START_PATTERNS.items():
        end_pat = STEP_PATTERNS[step]
        t_start = t_end = None
        for ts, msg in entries:
            if t_start is None and re.search(start_pat, msg, re.IGNORECASE):
                t_start = ts
            if t_start and re.search(end_pat, msg, re.IGNORECASE):
                t_end = ts
                break
        if t_start and t_end:
            durations[step] = (t_end - t_start).total_seconds() / 60  # minutes
    return durations


def collect_all(log_dir: Path) -> pd.DataFrame:
    """Parse all per-sample log files in log_dir."""
    records = []
    for log_file in sorted(log_dir.glob("*.log")):
        # Skip pipeline-wide and error logs
        if log_file.stem.startswith("analysis_") or log_file.stem.startswith("errors_"):
            continue
        # Sample name = everything before the last underscore+timestamp
        sample = re.sub(r"_\d{8}_\d{6}$", "", log_file.stem)
        durations = parse_log(log_file)
        if durations:
            row = {"sample": sample}
            row.update(durations)
            records.append(row)
    if not records:
        raise SystemExit(
            f"No per-sample log files found in {log_dir}.\n"
            "Run the pipeline first: exomeflow run ..."
        )
    return pd.DataFrame(records).set_index("sample")


def plot_heatmap(df: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(14, max(4, len(df) * 0.6)))
    sns.heatmap(
        df, annot=True, fmt=".1f", cmap="YlOrRd",
        linewidths=0.5, ax=ax,
        cbar_kws={"label": "Minutes"},
    )
    ax.set_title("ExomeFlow — Runtime per Step per Sample (minutes)", fontsize=13, pad=12)
    ax.set_xlabel("Pipeline Step")
    ax.set_ylabel("Sample")
    plt.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"  Saved: {out}")


def plot_bars(df: pd.DataFrame, out: Path) -> None:
    totals = df.sum(axis=1).sort_values(ascending=True)
    fig, axes = plt.subplots(1, 2, figsize=(16, max(4, len(df) * 0.5)))

    # Left — total runtime per sample
    totals.plot.barh(ax=axes[0], color="#1f77b4", edgecolor="white")
    axes[0].set_xlabel("Total Runtime (minutes)")
    axes[0].set_title("Total Runtime per Sample")
    axes[0].xaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f min"))
    for bar, val in zip(axes[0].patches, totals):
        axes[0].text(val + 0.3, bar.get_y() + bar.get_height() / 2,
                     f"{val:.1f}", va="center", fontsize=8)

    # Right — stacked bar per step
    step_order = list(START_PATTERNS.keys())
    cols = [c for c in step_order if c in df.columns]
    df_sorted = df.loc[totals.index, cols]
    df_sorted.plot.barh(stacked=True, ax=axes[1], colormap="tab10", edgecolor="white")
    axes[1].set_xlabel("Runtime (minutes)")
    axes[1].set_title("Step Breakdown per Sample")
    axes[1].legend(loc="lower right", fontsize=7)

    fig.suptitle("ExomeFlow — Performance Benchmark", fontsize=14, y=1.01)
    plt.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description="ExomeFlow performance benchmark")
    ap.add_argument("--log-dir", required=True, type=Path,
                    help="ExomeFlow logs directory (e.g. results/logs)")
    ap.add_argument("--out-dir", default="benchmarks/results", type=Path,
                    help="Output directory for tables and figures")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nParsing logs from: {args.log_dir}")
    df = collect_all(args.log_dir)

    # Fill missing steps with 0
    all_steps = list(START_PATTERNS.keys())
    for step in all_steps:
        if step not in df.columns:
            df[step] = 0.0
    df = df[[s for s in all_steps if s in df.columns]]

    # Summary
    print(f"\n{'Sample':<20} {'Total (min)':>12}")
    print("-" * 34)
    for sample, row in df.iterrows():
        print(f"  {sample:<18} {row.sum():>10.1f}")
    print(f"\n  Mean total runtime : {df.sum(axis=1).mean():.1f} min")
    print(f"  Max  total runtime : {df.sum(axis=1).max():.1f} min")
    print(f"  Slowest step (mean): {df.mean().idxmax()} ({df.mean().max():.1f} min)\n")

    # Save CSV
    csv_path = args.out_dir / "runtime_table.csv"
    df.to_csv(csv_path, float_format="%.2f")
    print(f"  Saved: {csv_path}")

    # Figures
    plot_heatmap(df, args.out_dir / "runtime_heatmap.png")
    plot_bars(df, args.out_dir / "runtime_bars.png")
    print("\nDone.")


if __name__ == "__main__":
    main()
