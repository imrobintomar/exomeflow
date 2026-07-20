"""
Read-depth CNV calling — GATK CollectReadCounts -> DenoiseReadCounts ->
PlotDenoisedCopyRatios, single-sample, no panel-of-normals required.

Gated on cfg.call_cnv (see steps.py). Coarser than a proper cohort PoN-based
gCNV call, but usable without requiring users to build a panel of normals
first. Requires --intervals (read-depth binning needs bounded regions).
"""

from __future__ import annotations

import logging
import shutil
from typing import TYPE_CHECKING

from exomeflow.utils import Checkpoint, PipelineStepError, run_cmd

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")

STEP = "cnv"


def run_cnv_calling(sample: str, cfg: "Config", checkpoint: Checkpoint) -> None:
    """
    Call read-depth CNVs for *sample*.

    Input  : <map_dir>/<sample>_recalibrated.bam
    Output : <cnv_dir>/<sample>_denoised_cr.tsv + .png plot
    """
    if checkpoint.done(sample, STEP):
        logger.info("[%s] CNV calling already completed, skipping.", sample)
        return

    if not cfg.has_intervals:
        raise PipelineStepError(
            f"[{sample}] --cnv requires --intervals "
            "(read-depth CNV calling needs bounded genomic bins)."
        )

    env = cfg.env()
    bam = cfg.map_dir / f"{sample}_recalibrated.bam"
    counts = cfg.cnv_dir / f"{sample}_counts.hdf5"
    standardized = cfg.cnv_dir / f"{sample}_standardized_cr.tsv"
    denoised = cfg.cnv_dir / f"{sample}_denoised_cr.tsv"
    seq_dict = cfg.reference.with_suffix(".dict")

    logger.info("[%s] Running CollectReadCounts ...", sample)
    run_cmd(
        [
            "gatk", "CollectReadCounts",
            "-I", str(bam),
            "-L", str(cfg.intervals),
            "--interval-merging-rule", "OVERLAPPING_ONLY",
            "-O", str(counts),
            "--format", "HDF5",
        ],
        env=env, step_name="CollectReadCounts", sample=sample,
    )

    logger.info("[%s] Running DenoiseReadCounts (no panel-of-normals) ...", sample)
    run_cmd(
        [
            "gatk", "DenoiseReadCounts",
            "-I", str(counts),
            "--standardized-copy-ratios", str(standardized),
            "--denoised-copy-ratios", str(denoised),
        ],
        env=env, step_name="DenoiseReadCounts", sample=sample,
    )

    # GATK's PlotDenoisedCopyRatios shells out to Rscript (optparse/
    # data.table/ggplot2) and, per GATK's own documented behavior, degrades
    # gracefully (exit 0, no plot produced) rather than hard-failing when R
    # isn't available — so its absence would otherwise be invisible here.
    # Found via audit: no R dependency check exists anywhere in this codebase.
    if not shutil.which("Rscript"):
        logger.warning(
            "[%s] Rscript not found on PATH — PlotDenoisedCopyRatios may run "
            "but produce no plot (install R with optparse/data.table/ggplot2 "
            "to get the copy-ratio plot).",
            sample,
        )

    logger.info("[%s] Running PlotDenoisedCopyRatios ...", sample)
    run_cmd(
        [
            "gatk", "PlotDenoisedCopyRatios",
            "--standardized-copy-ratios", str(standardized),
            "--denoised-copy-ratios", str(denoised),
            "--sequence-dictionary", str(seq_dict),
            "--output", str(cfg.cnv_dir),
            "--output-prefix", sample,
        ],
        env=env, step_name="PlotDenoisedCopyRatios", sample=sample,
    )

    # A reported-successful exit doesn't guarantee real output (same class
    # of bug already fixed for SortSam/ApplyBQSR/ANNOVAR elsewhere in this
    # codebase) — the core deliverable (denoised copy ratios) must exist
    # before this step is checkpointed done, even though the plot itself may
    # legitimately be absent if Rscript is missing. Found via audit.
    if not (denoised.exists() and denoised.stat().st_size > 0):
        raise PipelineStepError(
            f"[{sample}] CNV calling reported success but {denoised} is missing or empty."
        )

    checkpoint.mark(sample, STEP)
    logger.log(25, "[%s] CNV calling completed: %s", sample, denoised)
