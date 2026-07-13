"""
Cohort step — joint genotyping across all samples' per-sample GVCFs.

Only runs when `cfg.joint_genotyping` is True (opt-in; see steps.py). Mirrors
GATK's cohort best-practice: GenomicsDBImport -> GenotypeGVCFs, replacing the
per-sample VCF that `variant_calling.run_haplotype_caller` would otherwise
produce with one shared multi-sample GVCF workspace.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from exomeflow.utils import PipelineStepError, run_cmd

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")

STEP = "joint_genotyping"


def run_joint_genotyping(samples: list[str], cfg: "Config") -> None:
    """
    Combine per-sample GVCFs into one cohort VCF.

    Input  : <vcf_dir>/<sample>.g.vcf.gz for every sample
    Output : <cohort_dir>/cohort.vcf.gz
    """
    if not cfg.has_intervals:
        raise PipelineStepError(
            "[cohort] --joint-genotyping requires --intervals "
            "(GenomicsDBImport needs bounded genomic regions)."
        )

    env = cfg.env()
    workspace = cfg.cohort_dir / "genomicsdb"
    cohort_vcf = cfg.cohort_dir / "cohort.vcf.gz"

    logger.info("[cohort] Importing %d sample GVCF(s) into GenomicsDB ...", len(samples))

    import_cmd = ["gatk", "GenomicsDBImport"]
    for sample in samples:
        import_cmd += ["-V", str(cfg.vcf_dir / f"{sample}.g.vcf.gz")]
    import_cmd += [
        "--genomicsdb-workspace-path", str(workspace),
        "--overwrite-existing-genomicsdb-workspace", "true",
        "-L", str(cfg.intervals),
        "--interval-padding", str(cfg.interval_padding),
    ]
    run_cmd(import_cmd, env=env, step_name="GenomicsDBImport", sample="cohort")

    logger.info("[cohort] Running GenotypeGVCFs ...")

    run_cmd(
        [
            "gatk", "GenotypeGVCFs",
            "-R", str(cfg.reference),
            "-V", f"gendb://{workspace}",
            "-O", str(cohort_vcf),
            "--dbsnp", str(cfg.dbsnp),
        ],
        env=env, step_name="GenotypeGVCFs", sample="cohort",
    )

    logger.log(25, "[cohort] Joint genotyping completed: %s", cohort_vcf)
