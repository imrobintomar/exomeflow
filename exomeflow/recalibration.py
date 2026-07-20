"""
Step 7 — Base Quality Score Recalibration (BQSR).

Mirrors the Bash ``run_base_recalibration`` function:
  1. GATK BaseRecalibrator  → recal table
  2. GATK ApplyBQSR         → recalibrated BAM
  3. samtools index         → BAI
  4. Remove intermediate files (sorted, markdup, recal table)
     Keeps _recalibrated.bam + .bai for IGV inspection.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from exomeflow.utils import Checkpoint, PipelineStepError, run_cmd

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")

STEP = "bqsr"


def run_base_recalibration(
    sample: str, cfg: "Config", checkpoint: Checkpoint
) -> None:
    """
    Run BaseRecalibrator + ApplyBQSR for *sample*.

    Input  : <map_dir>/<sample>_markdup.bam
    Output : <map_dir>/<sample>_recalibrated.bam
             <map_dir>/<sample>_recalibrated.bam.bai
    Removes: _sorted.bam, _markdup.bam, _markdup.bai, recal.table
    """
    if checkpoint.done(sample, STEP):
        logger.info("[%s] BQSR already completed, skipping.", sample)
        return

    markdup_bam = cfg.map_dir / f"{sample}_markdup.bam"
    recal_table = cfg.map_dir / f"{sample}_recal.table"
    recal_bam = cfg.map_dir / f"{sample}_recalibrated.bam"

    # ------------------------------------------------------------------ 1/2
    logger.info("[%s] Running BaseRecalibrator ...", sample)

    run_cmd(
        [
            "gatk", "BaseRecalibrator",
            "-I", str(markdup_bam),
            "--known-sites", str(cfg.dbsnp),
            "--known-sites", str(cfg.known_indels),
            "--known-sites", str(cfg.mills),
            "-O", str(recal_table),
            "-R", str(cfg.reference),
        ],
        env=cfg.env(),
        step_name="BaseRecalibrator",
        sample=sample,
    )

    # ------------------------------------------------------------------ 2/2
    logger.info("[%s] Applying BQSR ...", sample)

    run_cmd(
        [
            "gatk", "ApplyBQSR",
            "--bqsr-recal-file", str(recal_table),
            "-I",                str(markdup_bam),
            "-O",                str(recal_bam),
            "-R",                str(cfg.reference),
        ],
        env=cfg.env(),
        step_name="ApplyBQSR",
        sample=sample,
    )

    # ------------------------------------------------------------------ index
    run_cmd(
        ["samtools", "index", str(recal_bam)],
        env=cfg.env(),
        step_name="samtools index (recalibrated BAM)",
        sample=sample,
    )

    # A reported-successful exit doesn't guarantee real output (same class
    # of bug already fixed for ANNOVAR/gsutil elsewhere in this codebase) —
    # verify before deleting the last remaining upstream copies. Found via
    # audit: a truncated/empty recalibrated BAM used to destroy sorted.bam/
    # markdup.bam/the recal table all at once, forcing a full re-align
    # instead of just retrying BQSR.
    recal_bai = cfg.map_dir / f"{sample}_recalibrated.bam.bai"
    if not (
        recal_bam.exists() and recal_bam.stat().st_size > 0
        and recal_bai.exists() and recal_bai.stat().st_size > 0
    ):
        raise PipelineStepError(
            f"[{sample}] BQSR reported success but {recal_bam.name}/{recal_bai.name} "
            f"is missing or empty — not deleting upstream intermediates."
        )

    # ------------------------------------------------------------------ clean
    # Keep _recalibrated.bam + .bai for IGV variant validation
    for p in [
        cfg.map_dir / f"{sample}_sorted.bam",
        cfg.map_dir / f"{sample}_markdup.bam",
        cfg.map_dir / f"{sample}_markdup.bai",
        recal_table,
    ]:
        p.unlink(missing_ok=True)

    logger.info(
        "[%s] Kept for IGV: %s",
        sample,
        cfg.map_dir / f"{sample}_recalibrated.bam",
    )

    checkpoint.mark(sample, STEP)
    logger.success("[%s] BQSR completed.", sample)
