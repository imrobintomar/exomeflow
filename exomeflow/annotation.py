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

from exomeflow.utils import Checkpoint, PipelineStepError, count_variants, run_cmd

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")

STEP = "annovar"


def annotate(
    *, label: str, input_vcf: Path, output_prefix: Path, cfg: "Config"
) -> None:
    """
    Build-agnostic ANNOVAR invocation shared by per-sample and cohort
    annotation. Produces `<output_prefix>.<buildver>_multianno.{vcf,txt}`.

    Skips cleanly (rather than crashing) when *input_vcf* has zero variants —
    ANNOVAR's table_annovar.pl doesn't handle a fully empty query gracefully
    (crashes with "the last column in header row should start with
    'Otherinfo'"), which a real sample can hit on catastrophic sequencing
    failure or an intervals BED with near-zero overlap with actual coverage.
    """
    if count_variants(input_vcf) == 0:
        logger.warning(
            "[%s] %s has 0 variants — skipping ANNOVAR annotation "
            "(nothing to annotate).",
            label, input_vcf,
        )
        return

    table_annovar = Path(cfg.annovar_bin) / "table_annovar.pl"

    logger.info("[%s] Running ANNOVAR annotation ...", label)

    cmd = [
        "perl", str(table_annovar),
        str(input_vcf),
        str(cfg.annovar_db),
        "--buildver", cfg.annovar_buildver,
        "--out",      str(output_prefix),
        "--remove",
        "--protocol", cfg.annovar_protocols,
        "--operation", cfg.annovar_operations,
        "-nastring", ".",
        "--polish",
        "--otherinfo",
        "--vcfinput",
        "--thread", str(cfg.annovar_threads),
    ]

    run_cmd(cmd, env=cfg.env(), step_name="table_annovar.pl", sample=label)

    # table_annovar.pl exiting 0 doesn't guarantee it actually wrote the
    # output — verified elsewhere in this codebase that tool exit codes
    # alone aren't trustworthy. A real (non-zero-variant) run must produce
    # both files; treat a missing one as a genuine failure rather than
    # silently letting the caller checkpoint a step with no real output.
    out_txt = Path(f"{output_prefix}.{cfg.annovar_buildver}_multianno.txt")
    out_vcf = Path(f"{output_prefix}.{cfg.annovar_buildver}_multianno.vcf")
    if not (out_txt.exists() and out_vcf.exists()):
        raise PipelineStepError(
            f"[{label}] table_annovar.pl exited 0 but expected output "
            f"({out_txt.name} / {out_vcf.name}) wasn't produced."
        )

    # Remove ANNOVAR intermediate file (--remove doesn't always clean this up)
    avinput = Path(f"{output_prefix}.avinput")
    avinput.unlink(missing_ok=True)


def run_annovar_annotation(
    sample: str, cfg: "Config", checkpoint: Checkpoint
) -> None:
    """
    Annotate PASS variants for *sample* using ANNOVAR table_annovar.pl.

    Input  : <vcf_dir>/<sample>_PASS.vcf
    Output : <vcf_dir>/<sample>.annovar.<buildver>_multianno.vcf
             <vcf_dir>/<sample>.annovar.<buildver>_multianno.txt
    """
    if checkpoint.done(sample, STEP):
        logger.info("[%s] ANNOVAR annotation already completed, skipping.", sample)
        return

    annotate(
        label=sample,
        input_vcf=cfg.vcf_dir / f"{sample}_PASS.vcf",
        output_prefix=cfg.vcf_dir / f"{sample}.annovar",
        cfg=cfg,
    )

    checkpoint.mark(sample, STEP)
    logger.log(25, "[%s] ANNOVAR annotation completed.", sample)


def run_cohort_annotation(samples: list[str], cfg: "Config") -> None:
    """
    Cohort-level counterpart of `run_annovar_annotation`, applied once to the
    joint-genotyped cohort PASS VCF instead of once per sample.

    Input  : <cohort_dir>/cohort_PASS.vcf
    Output : <cohort_dir>/cohort.annovar.<buildver>_multianno.{vcf,txt}
    """
    annotate(
        label="cohort",
        input_vcf=cfg.cohort_dir / "cohort_PASS.vcf",
        output_prefix=cfg.cohort_dir / "cohort.annovar",
        cfg=cfg,
    )
    logger.log(25, "[cohort] ANNOVAR annotation completed.")
