"""
Step 12 — Variant annotation with ANNOVAR.

Mirrors the Bash ``run_annovar_annotation`` function exactly:
  - Calls table_annovar.pl with protocols / operations from Config
  - Removes the intermediate .avinput file
  - Passes --thread for parallel ANNOVAR processing
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

from exomeflow.utils import Checkpoint, run_cmd

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")

STEP = "annovar"


def run_annovar_annotation(
    sample: str, cfg: "Config", checkpoint: Checkpoint
) -> None:
    """
    Annotate PASS variants for *sample* using ANNOVAR table_annovar.pl.

    Input  : <vcf_dir>/<sample>_PASS.vcf
    Output : <vcf_dir>/<sample>.annovar.hg38_multianno.vcf
             <vcf_dir>/<sample>.annovar.hg38_multianno.txt
    Removes: <vcf_dir>/<sample>.annovar.avinput (ANNOVAR intermediate)
    """
    if checkpoint.done(sample, STEP):
        logger.info("[%s] ANNOVAR annotation already completed, skipping.", sample)
        return

    vcf = cfg.vcf_dir / f"{sample}_PASS.vcf"
    prefix = cfg.vcf_dir / f"{sample}.annovar"

    table_annovar = Path(cfg.annovar_bin) / "table_annovar.pl"

    logger.info("[%s] Running ANNOVAR annotation ...", sample)

    cmd = [
        "perl", str(table_annovar),
        str(vcf),
        str(cfg.annovar_db),
        "--buildver", "hg38",
        "--out",      str(prefix),
        "--remove",
        "--protocol", cfg.annovar_protocols,
        "--operation", cfg.annovar_operations,
        "-nastring", ".",
        "--polish",
        "--otherinfo",
        "--vcfinput",
        "--thread", str(cfg.annovar_threads),
    ]

    run_cmd(cmd, env=cfg.env(), step_name="table_annovar.pl", sample=sample)

    # Remove ANNOVAR intermediate file (--remove doesn't always clean this up)
    avinput = cfg.vcf_dir / f"{sample}.annovar.avinput"
    avinput.unlink(missing_ok=True)

    checkpoint.mark(sample, STEP)
    logger.log(25, "[%s] ANNOVAR annotation completed.", sample)
