"""
Somatic mode (tumor-only) — Mutect2 calling + Mutect2's own filtering chain.

Replaces variant_calling.run_haplotype_caller / filtering.run_variant_filtration
when cfg.mode == "somatic" (mutually exclusive via steps.py applies() gates).
Tumor-normal pairing is out of scope for V2 — every detected sample is called
tumor-only; --germline-resource (gnomAD AF-only VCF) is optional but strongly
recommended to keep the false-positive rate down without a matched normal.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from exomeflow.utils import Checkpoint, run_cmd

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")

MUTECT2_STEP = "mutect2"
FILTER_STEP = "somatic_filter"


def run_mutect2(sample: str, cfg: "Config", checkpoint: Checkpoint) -> None:
    """
    Call somatic variants tumor-only for *sample* with Mutect2.

    Input  : <map_dir>/<sample>_recalibrated.bam
    Output : <vcf_dir>/<sample>_unfiltered.vcf.gz
    """
    if checkpoint.done(sample, MUTECT2_STEP):
        logger.info("[%s] Mutect2 already completed, skipping.", sample)
        return

    bam = cfg.map_dir / f"{sample}_recalibrated.bam"
    unfiltered = cfg.vcf_dir / f"{sample}_unfiltered.vcf.gz"

    if not cfg.germline_resource:
        logger.warning(
            "[%s] No --germline-resource supplied — tumor-only Mutect2 calls "
            "will have a higher false-positive rate without a population AF prior.",
            sample,
        )

    cmd = [
        "gatk", "Mutect2",
        "-R", str(cfg.reference),
        "-I", str(bam),
        "-tumor", sample,
        "-O", str(unfiltered),
    ]
    if cfg.germline_resource:
        cmd += ["--germline-resource", str(cfg.germline_resource)]
    if cfg.has_intervals:
        cmd += ["-L", str(cfg.intervals), "--interval-padding", str(cfg.interval_padding)]

    logger.info("[%s] Running Mutect2 (tumor-only) ...", sample)
    run_cmd(cmd, env=cfg.env(), step_name="Mutect2", sample=sample)

    checkpoint.mark(sample, MUTECT2_STEP)
    logger.log(25, "[%s] Mutect2 completed.", sample)


def run_somatic_filtration(sample: str, cfg: "Config", checkpoint: Checkpoint) -> None:
    """
    Estimate contamination and apply Mutect2's own filtering chain.

    Input  : <vcf_dir>/<sample>_unfiltered.vcf.gz
    Output : <vcf_dir>/<sample>_PASS.vcf
    """
    if checkpoint.done(sample, FILTER_STEP):
        logger.info("[%s] Somatic filtering already completed, skipping.", sample)
        return

    env = cfg.env()
    bam = cfg.map_dir / f"{sample}_recalibrated.bam"
    unfiltered = cfg.vcf_dir / f"{sample}_unfiltered.vcf.gz"
    filtered = cfg.vcf_dir / f"{sample}_filtered.vcf.gz"
    pass_vcf = cfg.vcf_dir / f"{sample}_PASS.vcf"

    filter_cmd = [
        "gatk", "FilterMutectCalls",
        "-R", str(cfg.reference),
        "-V", str(unfiltered),
        "-O", str(filtered),
    ]

    if cfg.germline_resource:
        pileups = cfg.vcf_dir / f"{sample}_pileups.table"
        contamination = cfg.vcf_dir / f"{sample}_contamination.table"

        logger.info("[%s] Running GetPileupSummaries ...", sample)
        pileup_cmd = [
            "gatk", "GetPileupSummaries",
            "-I", str(bam),
            "-V", str(cfg.germline_resource),
            "-L", str(cfg.germline_resource),
            "-O", str(pileups),
        ]
        run_cmd(pileup_cmd, env=env, step_name="GetPileupSummaries", sample=sample)

        logger.info("[%s] Running CalculateContamination ...", sample)
        run_cmd(
            ["gatk", "CalculateContamination",
             "-I", str(pileups), "-O", str(contamination)],
            env=env, step_name="CalculateContamination", sample=sample,
        )
        filter_cmd += ["--contamination-table", str(contamination)]

    logger.info("[%s] Running FilterMutectCalls ...", sample)
    run_cmd(filter_cmd, env=env, step_name="FilterMutectCalls", sample=sample)

    logger.info("[%s] Extracting PASS variants ...", sample)
    run_cmd(
        ["gatk", "SelectVariants",
         "-R", str(cfg.reference),
         "-V", str(filtered),
         "-O", str(pass_vcf),
         "--exclude-filtered",
         "--exclude-non-variants"],
        env=env, step_name="SelectVariants (PASS)", sample=sample,
    )

    checkpoint.mark(sample, FILTER_STEP)
    logger.log(25, "[%s] Somatic filtering completed.", sample)
