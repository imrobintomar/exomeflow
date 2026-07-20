"""
Steps 9-11 — Variant filtering (GATK hard filters).

Mirrors the Bash ``run_variant_filtration`` function exactly:
  1. SelectVariants  → separate SNPs
  2. SelectVariants  → separate INDELs
  3. VariantFiltration (SNP thresholds)
  4. VariantFiltration (INDEL thresholds)
  5. MergeVcfs       → merge filtered SNPs + INDELs
  6. SelectVariants  → PASS-only extraction (--exclude-filtered --exclude-non-variants)
  7. Summary log (total / SNP / INDEL / PASS counts)
  8. Remove intermediate VCFs (keep raw + PASS)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

from exomeflow.config import (
    INDEL_FILTERS,
    INDEL_GENOTYPE_FILTERS,
    SNP_FILTERS,
    SNP_GENOTYPE_FILTERS,
)
from exomeflow.utils import Checkpoint, PipelineStepError, count_variants, run_cmd

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")

STEP = "filter"


def hard_filter(
    *,
    label: str,
    input_vcf: Path,
    output_pass_vcf: Path,
    workdir: Path,
    reference: Path,
    env: dict[str, str],
) -> None:
    """
    Build-agnostic GATK hard-filter chain shared by the per-sample and
    cohort filtration paths (`run_variant_filtration` / cohort_filtering.py).

    Input  : *input_vcf* (either a single-sample or cohort-joint-genotyped VCF)
    Output : *output_pass_vcf*
    """
    snp_raw = workdir / f"{label}_snp_raw.vcf"
    indel_raw = workdir / f"{label}_indel_raw.vcf"
    snp_filtered = workdir / f"{label}_snp_filtered.vcf"
    indel_filtered = workdir / f"{label}_indel_filtered.vcf"
    merged = workdir / f"{label}_merged_filtered.vcf"

    # ------------------------------------------------------------------ 1. Separate
    logger.info("[%s] Separating SNPs and INDELs ...", label)

    run_cmd(
        ["gatk", "SelectVariants",
         "-R", str(reference), "-V", str(input_vcf),
         "--select-type-to-include", "SNP",
         "-O", str(snp_raw)],
        env=env, step_name="SelectVariants (SNP)", sample=label,
    )

    run_cmd(
        ["gatk", "SelectVariants",
         "-R", str(reference), "-V", str(input_vcf),
         "--select-type-to-include", "INDEL",
         "-O", str(indel_raw)],
        env=env, step_name="SelectVariants (INDEL)", sample=label,
    )

    # ------------------------------------------------------------------ 2. Filter SNPs
    logger.info("[%s] Filtering SNPs ...", label)

    snp_cmd = [
        "gatk", "VariantFiltration",
        "-R", str(reference),
        "-V", str(snp_raw),
        "-O", str(snp_filtered),
    ]
    for expr, name in SNP_FILTERS:
        snp_cmd += ["--filter-expression", expr, "--filter-name", name]
    for expr, name in SNP_GENOTYPE_FILTERS:
        snp_cmd += ["--genotype-filter-expression", expr,
                    "--genotype-filter-name", name]

    run_cmd(snp_cmd, env=env, step_name="VariantFiltration (SNP)", sample=label)

    # ------------------------------------------------------------------ 3. Filter INDELs
    logger.info("[%s] Filtering INDELs ...", label)

    indel_cmd = [
        "gatk", "VariantFiltration",
        "-R", str(reference),
        "-V", str(indel_raw),
        "-O", str(indel_filtered),
    ]
    for expr, name in INDEL_FILTERS:
        indel_cmd += ["--filter-expression", expr, "--filter-name", name]
    for expr, name in INDEL_GENOTYPE_FILTERS:
        indel_cmd += ["--genotype-filter-expression", expr,
                      "--genotype-filter-name", name]

    run_cmd(indel_cmd, env=env, step_name="VariantFiltration (INDEL)", sample=label)

    # ------------------------------------------------------------------ 4. Merge
    logger.info("[%s] Merging filtered SNPs and INDELs ...", label)

    run_cmd(
        ["gatk", "MergeVcfs",
         "-I", str(snp_filtered),
         "-I", str(indel_filtered),
         "-O", str(merged)],
        env=env, step_name="MergeVcfs", sample=label,
    )

    # ------------------------------------------------------------------ 5. PASS
    logger.info("[%s] Extracting PASS variants ...", label)

    run_cmd(
        ["gatk", "SelectVariants",
         "-R", str(reference),
         "-V", str(merged),
         "-O", str(output_pass_vcf),
         "--exclude-filtered",
         "--exclude-non-variants"],
        env=env, step_name="SelectVariants (PASS)", sample=label,
    )

    # A reported-successful exit doesn't guarantee real output (same class
    # of bug already fixed for SortSam/ApplyBQSR/ANNOVAR elsewhere in this
    # codebase) — verify before deleting the intermediates below, and
    # before the caller checkpoints this step done. Found via audit.
    if not (output_pass_vcf.exists() and output_pass_vcf.stat().st_size > 0):
        raise PipelineStepError(
            f"[{label}] SelectVariants (PASS) reported success but "
            f"{output_pass_vcf} is missing or empty — not deleting intermediates."
        )

    # ------------------------------------------------------------------ 6. Summary
    total = count_variants(input_vcf)
    snp_count = count_variants(snp_raw)
    indel_count = count_variants(indel_raw)
    passed = count_variants(output_pass_vcf)
    logger.info(
        "[%s] Total: %d | SNPs: %d | INDELs: %d | PASS: %d | Filtered: %d",
        label, total, snp_count, indel_count, passed, total - passed,
    )

    # ------------------------------------------------------------------ 7. Clean-up
    intermediates = [
        snp_raw, Path(str(snp_raw) + ".idx"),
        indel_raw, Path(str(indel_raw) + ".idx"),
        snp_filtered, Path(str(snp_filtered) + ".idx"),
        indel_filtered, Path(str(indel_filtered) + ".idx"),
        merged, Path(str(merged) + ".idx"),
    ]
    for p in intermediates:
        p.unlink(missing_ok=True)


def run_variant_filtration(
    sample: str, cfg: "Config", checkpoint: Checkpoint
) -> None:
    """
    Apply GATK hard filters and extract PASS variants for *sample*.

    Input  : <vcf_dir>/<sample>.vcf
    Output : <vcf_dir>/<sample>_PASS.vcf
    Keeps  : <vcf_dir>/<sample>.vcf (raw HaplotypeCaller output)
    """
    if checkpoint.done(sample, STEP):
        logger.info("[%s] Variant filtering already completed, skipping.", sample)
        return

    hard_filter(
        label=sample,
        input_vcf=cfg.vcf_dir / f"{sample}.vcf",
        output_pass_vcf=cfg.vcf_dir / f"{sample}_PASS.vcf",
        workdir=cfg.vcf_dir,
        reference=cfg.reference,
        env=cfg.env(),
    )

    checkpoint.mark(sample, STEP)
    logger.log(25, "[%s] Variant filtering completed.", sample)


def run_cohort_filtration(samples: list[str], cfg: "Config") -> None:
    """
    Cohort-level counterpart of `run_variant_filtration`, applied once to the
    joint-genotyped cohort VCF instead of once per sample.

    Input  : <cohort_dir>/cohort.vcf.gz
    Output : <cohort_dir>/cohort_PASS.vcf
    """
    hard_filter(
        label="cohort",
        input_vcf=cfg.cohort_dir / "cohort.vcf.gz",
        output_pass_vcf=cfg.cohort_dir / "cohort_PASS.vcf",
        workdir=cfg.cohort_dir,
        reference=cfg.reference,
        env=cfg.env(),
    )
    logger.log(25, "[cohort] Variant filtering completed.")
