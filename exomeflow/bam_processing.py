"""
Steps 3-6 — BAM processing.

Mirrors the Bash functions:
  sort_bam           → GATK SortSam (coordinate order)
  generate_flagstat  → samtools flagstat
  mark_duplicates    → GATK MarkDuplicates
  build_bam_index    → GATK BuildBamIndex
"""

from __future__ import annotations

import logging
import subprocess
from typing import TYPE_CHECKING

from exomeflow.utils import Checkpoint, PipelineStepError, run_cmd

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")


# ---------------------------------------------------------------------------
# Step 3 — SortSam
# ---------------------------------------------------------------------------

def sort_bam(sample: str, cfg: "Config", checkpoint: Checkpoint) -> None:
    """
    Sort *sample*.bam by coordinate using GATK SortSam.

    Input  : <map_dir>/<sample>.bam
    Output : <map_dir>/<sample>_sorted.bam
    Removes: <map_dir>/<sample>.bam  (raw unsorted BAM freed after sorting)
    """
    if checkpoint.done(sample, "sort"):
        logger.info("[%s] BAM sorting already completed, skipping.", sample)
        return

    input_bam = cfg.map_dir / f"{sample}.bam"
    output_bam = cfg.map_dir / f"{sample}_sorted.bam"

    logger.info("[%s] Sorting BAM ...", sample)

    run_cmd(
        ["gatk", "SortSam",
         "-I", str(input_bam),
         "-O", str(output_bam),
         "-SO", "coordinate"],
        env=cfg.env(),
        step_name="SortSam",
        sample=sample,
    )

    # A reported-successful exit doesn't guarantee real output (same class
    # of bug already fixed for ANNOVAR/gsutil elsewhere in this codebase) —
    # verify before deleting the *only* upstream copy. Found via audit: a
    # truncated/empty sorted BAM from a disk-full or killed-mid-write JVM
    # used to destroy the sole recovery path, forcing a full BWA re-align
    # instead of just retrying this one step.
    if not (output_bam.exists() and output_bam.stat().st_size > 0):
        raise PipelineStepError(
            f"[{sample}] SortSam exited 0 but {output_bam} is missing or empty "
            f"— not deleting {input_bam}."
        )

    input_bam.unlink(missing_ok=True)

    checkpoint.mark(sample, "sort")
    logger.log(25, "[%s] BAM sorting completed.", sample)


# ---------------------------------------------------------------------------
# Step 4 — flagstat
# ---------------------------------------------------------------------------

def generate_flagstat(sample: str, cfg: "Config", checkpoint: Checkpoint) -> None:
    """
    Run samtools flagstat on the sorted BAM.

    Input  : <map_dir>/<sample>_sorted.bam
    Output : <map_dir>/<sample>_flagstat.txt
    """
    if checkpoint.done(sample, "flagstat"):
        logger.info("[%s] Flagstat already generated, skipping.", sample)
        return

    bam = cfg.map_dir / f"{sample}_sorted.bam"
    output = cfg.map_dir / f"{sample}_flagstat.txt"

    logger.info("[%s] Generating flagstat ...", sample)

    # samtools flagstat writes to stdout; capture and write to file
    env = cfg.env()
    result = subprocess.run(
        ["samtools", "flagstat", str(bam)],
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    if result.returncode != 0:
        from exomeflow.utils import PipelineStepError
        raise PipelineStepError(
            f"[{sample}] samtools flagstat failed:\n{result.stderr}"
        )

    output.write_text(result.stdout, encoding="utf-8")

    checkpoint.mark(sample, "flagstat")
    logger.log(25, "[%s] Flagstat completed.", sample)


# ---------------------------------------------------------------------------
# Step 5 — MarkDuplicates
# ---------------------------------------------------------------------------

def mark_duplicates(sample: str, cfg: "Config", checkpoint: Checkpoint) -> None:
    """
    Mark PCR duplicates with GATK MarkDuplicates.

    Input  : <map_dir>/<sample>_sorted.bam
    Output : <map_dir>/<sample>_markdup.bam
             <map_dir>/<sample>_markdup_metrics.txt
    """
    if checkpoint.done(sample, "markdup"):
        logger.info("[%s] MarkDuplicates already completed, skipping.", sample)
        return

    input_bam = cfg.map_dir / f"{sample}_sorted.bam"
    output_bam = cfg.map_dir / f"{sample}_markdup.bam"
    metrics = cfg.map_dir / f"{sample}_markdup_metrics.txt"

    logger.info("[%s] Marking duplicates ...", sample)

    run_cmd(
        ["gatk", "MarkDuplicates",
         "-I", str(input_bam),
         "-O", str(output_bam),
         "-M", str(metrics)],
        env=cfg.env(),
        step_name="MarkDuplicates",
        sample=sample,
    )

    checkpoint.mark(sample, "markdup")
    logger.log(25, "[%s] MarkDuplicates completed.", sample)


# ---------------------------------------------------------------------------
# Step 6 — BuildBamIndex
# ---------------------------------------------------------------------------

def build_bam_index(sample: str, cfg: "Config", checkpoint: Checkpoint) -> None:
    """
    Build a BAI index for the duplicate-marked BAM.

    Input  : <map_dir>/<sample>_markdup.bam
    Output : <map_dir>/<sample>_markdup.bai
    """
    if checkpoint.done(sample, "index"):
        logger.info("[%s] BAM index already created, skipping.", sample)
        return

    bam = cfg.map_dir / f"{sample}_markdup.bam"

    logger.info("[%s] Building BAM index ...", sample)

    run_cmd(
        ["gatk", "BuildBamIndex", "-I", str(bam)],
        env=cfg.env(),
        step_name="BuildBamIndex",
        sample=sample,
    )

    checkpoint.mark(sample, "index")
    logger.log(25, "[%s] BAM index created.", sample)
