"""
Step 1 — Quality control with fastp.

Mirrors the Bash ``run_fastp`` function exactly:
  - min read length: 50 bp
  - min base quality: Q30
  - produces HTML + JSON QC reports
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from exomeflow.utils import Checkpoint, resolve_fastq_pair, run_cmd

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")

STEP = "fastp"


def run_fastp(sample: str, cfg: "Config", checkpoint: Checkpoint) -> None:
    """
    Run fastp quality-filtering for *sample*.

    Input
    -----
    <input_dir>/<sample>_1.fastq.gz  or  <input_dir>/<sample>_R1.fastq.gz
    <input_dir>/<sample>_2.fastq.gz  or  <input_dir>/<sample>_R2.fastq.gz

    Output
    ------
    <fastp_dir>/<sample>_1_filtered.fastq.gz
    <fastp_dir>/<sample>_2_filtered.fastq.gz
    <fastp_dir>/<sample>_fastp.html
    <fastp_dir>/<sample>_fastp.json
    """
    if checkpoint.done(sample, STEP):
        logger.info("[%s] fastp already completed, skipping.", sample)
        return

    r1, r2 = resolve_fastq_pair(cfg.input_dir, sample)

    out_r1 = cfg.fastp_dir / f"{sample}_1_filtered.fastq.gz"
    out_r2 = cfg.fastp_dir / f"{sample}_2_filtered.fastq.gz"
    html = cfg.fastp_dir / f"{sample}_fastp.html"
    json = cfg.fastp_dir / f"{sample}_fastp.json"

    logger.info("[%s] Running fastp ...", sample)

    cmd = [
        "fastp",
        "-i", str(r1),
        "-I", str(r2),
        "--length_required", "50",
        "--qualified_quality_phred", "30",
        "-o", str(out_r1),
        "-O", str(out_r2),
        "--html", str(html),
        "--json", str(json),
        "-w", str(cfg.fastp_threads),
    ]

    run_cmd(cmd, env=cfg.env(), step_name="fastp", sample=sample)

    checkpoint.mark(sample, STEP)
    logger.success("[%s] fastp completed.", sample)   # 25 = SUCCESS
