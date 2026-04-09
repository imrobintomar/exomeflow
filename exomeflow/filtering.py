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
from exomeflow.utils import Checkpoint, count_variants, run_cmd

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")

STEP = "filter"


def run_variant_filtration(
    sample: str, cfg: "Config", checkpoint: Checkpoint
) -> None:
    """
    Apply GATK hard filters and extract PASS variants for *sample*.

    Input  : <vcf_dir>/<sample>.vcf
    Output : <vcf_dir>/<sample>_PASS.vcf
    Keeps  : <vcf_dir>/<sample>.vcf (raw HaplotypeCaller output)
    Removes: snp_raw, indel_raw, snp_filtered, indel_filtered, merged_filtered
    """
    if checkpoint.done(sample, STEP):
        logger.info("[%s] Variant filtering already completed, skipping.", sample)
        return

    vcf           = cfg.vcf_dir / f"{sample}.vcf"
    snp_raw       = cfg.vcf_dir / f"{sample}_snp_raw.vcf"
    indel_raw     = cfg.vcf_dir / f"{sample}_indel_raw.vcf"
    snp_filtered  = cfg.vcf_dir / f"{sample}_snp_filtered.vcf"
    indel_filtered= cfg.vcf_dir / f"{sample}_indel_filtered.vcf"
    merged        = cfg.vcf_dir / f"{sample}_merged_filtered.vcf"
    pass_vcf      = cfg.vcf_dir / f"{sample}_PASS.vcf"

    env = cfg.env()

    # ------------------------------------------------------------------ 1. Separate
    logger.info("[%s] Separating SNPs and INDELs ...", sample)

    run_cmd(
        ["gatk", "SelectVariants",
         "-R", str(cfg.reference), "-V", str(vcf),
         "--select-type-to-include", "SNP",
         "-O", str(snp_raw)],
        env=env, step_name="SelectVariants (SNP)", sample=sample,
    )

    run_cmd(
        ["gatk", "SelectVariants",
         "-R", str(cfg.reference), "-V", str(vcf),
         "--select-type-to-include", "INDEL",
         "-O", str(indel_raw)],
        env=env, step_name="SelectVariants (INDEL)", sample=sample,
    )

    # ------------------------------------------------------------------ 2. Filter SNPs
    logger.info("[%s] Filtering SNPs ...", sample)

    snp_cmd = [
        "gatk", "VariantFiltration",
        "-R", str(cfg.reference),
        "-V", str(snp_raw),
        "-O", str(snp_filtered),
    ]
    for expr, name in SNP_FILTERS:
        snp_cmd += ["--filter-expression", expr, "--filter-name", name]
    for expr, name in SNP_GENOTYPE_FILTERS:
        snp_cmd += ["--genotype-filter-expression", expr,
                    "--genotype-filter-name", name]

    run_cmd(snp_cmd, env=env, step_name="VariantFiltration (SNP)", sample=sample)

    # ------------------------------------------------------------------ 3. Filter INDELs
    logger.info("[%s] Filtering INDELs ...", sample)

    indel_cmd = [
        "gatk", "VariantFiltration",
        "-R", str(cfg.reference),
        "-V", str(indel_raw),
        "-O", str(indel_filtered),
    ]
    for expr, name in INDEL_FILTERS:
        indel_cmd += ["--filter-expression", expr, "--filter-name", name]
    for expr, name in INDEL_GENOTYPE_FILTERS:
        indel_cmd += ["--genotype-filter-expression", expr,
                      "--genotype-filter-name", name]

    run_cmd(indel_cmd, env=env, step_name="VariantFiltration (INDEL)", sample=sample)

    # ------------------------------------------------------------------ 4. Merge
    logger.info("[%s] Merging filtered SNPs and INDELs ...", sample)

    run_cmd(
        ["gatk", "MergeVcfs",
         "-I", str(snp_filtered),
         "-I", str(indel_filtered),
         "-O", str(merged)],
        env=env, step_name="MergeVcfs", sample=sample,
    )

    # ------------------------------------------------------------------ 5. PASS
    logger.info("[%s] Extracting PASS variants ...", sample)

    run_cmd(
        ["gatk", "SelectVariants",
         "-R", str(cfg.reference),
         "-V", str(merged),
         "-O", str(pass_vcf),
         "--exclude-filtered",
         "--exclude-non-variants"],
        env=env, step_name="SelectVariants (PASS)", sample=sample,
    )

    # ------------------------------------------------------------------ 6. Summary
    total       = count_variants(vcf)
    snp_count   = count_variants(snp_raw)
    indel_count = count_variants(indel_raw)
    passed      = count_variants(pass_vcf)
    logger.info(
        "[%s] Total: %d | SNPs: %d | INDELs: %d | PASS: %d | Filtered: %d",
        sample, total, snp_count, indel_count, passed, total - passed,
    )

    # ------------------------------------------------------------------ 7. Clean-up
    intermediates = [
        snp_raw,    Path(str(snp_raw)    + ".idx"),
        indel_raw,  Path(str(indel_raw)  + ".idx"),
        snp_filtered, Path(str(snp_filtered)  + ".idx"),
        indel_filtered, Path(str(indel_filtered) + ".idx"),
        merged,     Path(str(merged)     + ".idx"),
    ]
    for p in intermediates:
        p.unlink(missing_ok=True)

    checkpoint.mark(sample, STEP)
    logger.log(25, "[%s] Variant filtering completed.", sample)
