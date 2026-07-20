"""
Step 2 — Read alignment with BWA MEM.

Mirrors the Bash ``run_bwa_mem`` function exactly:
  - BWA MEM flags: -Y -K 100000000
  - Read-group tag set to sample name (ID, PU, SM, LB = sample; PL = illumina)
  - Output piped through samtools view -Shb to produce a raw BAM
"""

from __future__ import annotations

import logging
import subprocess
from typing import TYPE_CHECKING

from exomeflow.utils import Checkpoint, PipelineStepError

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")

STEP = "bwa"


def run_bwa_mem(sample: str, cfg: "Config", checkpoint: Checkpoint) -> None:
    """
    Align filtered reads for *sample* with BWA MEM and convert to BAM.

    Input
    -----
    <fastp_dir>/<sample>_1_filtered.fastq.gz
    <fastp_dir>/<sample>_2_filtered.fastq.gz

    Output
    ------
    <map_dir>/<sample>.bam
    """
    if checkpoint.done(sample, STEP):
        logger.info("[%s] BWA MEM already completed, skipping.", sample)
        return

    r1 = cfg.fastp_dir / f"{sample}_1_filtered.fastq.gz"
    r2 = cfg.fastp_dir / f"{sample}_2_filtered.fastq.gz"
    output = cfg.map_dir / f"{sample}.bam"

    read_group = (
        f"@RG\\tID:{sample}\\tPU:{sample}\\tSM:{sample}"
        f"\\tLB:{sample}\\tPL:illumina"
    )

    logger.info("[%s] Running BWA MEM ...", sample)

    # BWA MEM → samtools view pipe (mirrors the Bash pipe)
    bwa_cmd = [
        "bwa", "mem",
        "-Y",
        "-K", "100000000",
        "-t", str(cfg.threads),
        "-R", read_group,
        str(cfg.reference),
        str(r1),
        str(r2),
    ]

    samtools_cmd = [
        "samtools", "view",
        "-Shb",
        "-o", str(output),
        "-",
    ]

    env = cfg.env()

    bwa_proc = subprocess.Popen(
        bwa_cmd,
        stdout=subprocess.PIPE,
        stderr=None,
        env=env,
    )

    try:
        samtools_proc = subprocess.Popen(
            samtools_cmd,
            stdin=bwa_proc.stdout,
            env=env,
        )
    except Exception:
        # If samtools fails to even start (e.g. missing from PATH), bwa_proc
        # is already spawned and would otherwise leak as an orphaned running
        # process with nothing left to read its output. Found via audit.
        bwa_proc.kill()
        bwa_proc.wait()
        raise

    # Allow bwa to receive SIGPIPE if samtools exits early
    if bwa_proc.stdout:
        bwa_proc.stdout.close()

    samtools_rc = samtools_proc.wait()
    bwa_rc = bwa_proc.wait()

    if bwa_rc != 0 or samtools_rc != 0:
        raise PipelineStepError(
            f"[{sample}] BWA MEM / samtools pipe failed "
            f"(bwa={bwa_rc}, samtools={samtools_rc})"
        )

    checkpoint.mark(sample, STEP)
    logger.success("[%s] BWA MEM and conversion to BAM completed.", sample)
