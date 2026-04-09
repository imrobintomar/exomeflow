"""
Step 8 — Variant calling with GATK HaplotypeCaller.

Mirrors the Bash ``run_haplotype_caller`` function:
  - Uses exome intervals BED + padding when provided; warns and falls back
    to whole-genome mode if the BED is absent.
  - native-pair-hmm-threads set to cfg.threads.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

from exomeflow.utils import Checkpoint, run_cmd

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")

STEP = "haplotype"


def run_haplotype_caller(
    sample: str, cfg: "Config", checkpoint: Checkpoint
) -> None:
    """
    Call variants with HaplotypeCaller for *sample*.

    Input  : <map_dir>/<sample>_recalibrated.bam
    Output : <vcf_dir>/<sample>.vcf
    """
    if checkpoint.done(sample, STEP):
        logger.info("[%s] HaplotypeCaller already completed, skipping.", sample)
        return

    bam = cfg.map_dir / f"{sample}_recalibrated.bam"
    vcf = cfg.vcf_dir / f"{sample}.vcf"

    logger.info("[%s] Running HaplotypeCaller ...", sample)

    cmd: list[str] = [
        "gatk", "HaplotypeCaller",
        "-R",     str(cfg.reference),
        "-I",     str(bam),
        "-O",     str(vcf),
        "--dbsnp", str(cfg.dbsnp),
        "--native-pair-hmm-threads", str(cfg.threads),
    ]

    # Exome-interval support
    if cfg.intervals and Path(cfg.intervals).exists():
        cmd += ["-L", str(cfg.intervals),
                "--interval-padding", str(cfg.interval_padding)]
        logger.info(
            "[%s] Using exome intervals: %s (padding: %d bp)",
            sample, cfg.intervals, cfg.interval_padding,
        )
    else:
        logger.warning(
            "[%s] No exome intervals BED found at '%s' — "
            "calling whole genome (slower, more false positives).",
            sample, cfg.intervals,
        )

    run_cmd(cmd, env=cfg.env(), step_name="HaplotypeCaller", sample=sample)

    checkpoint.mark(sample, STEP)
    logger.log(25, "[%s] HaplotypeCaller completed.", sample)
