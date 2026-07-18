"""
Orchestration layer — processes a single sample end-to-end and supports
parallel cohort processing via concurrent.futures.ProcessPoolExecutor.

SAMPLE_STEPS / COHORT_STEPS (steps.py) replace what used to be a hardcoded
step tuple: each step carries an `applies(cfg)` predicate so germline vs.
somatic vs. joint-genotyping vs. CNV all compose from the same list instead
of branching scattered through this file. See the V2 plan's "non-negotiable
invariant": per-sample annotated output is the default for any number of
samples; cohort (joint-genotyping) output only replaces it when the user
explicitly opts in via --joint-genotyping.
"""

from __future__ import annotations

import logging
import signal
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from typing import TYPE_CHECKING

from exomeflow.acmg_classification import run_intervar, run_intervar_cohort
from exomeflow.alignment import run_bwa_mem
from exomeflow.annotation import run_annovar_annotation, run_cohort_annotation
from exomeflow.bam_processing import (
    build_bam_index,
    generate_flagstat,
    mark_duplicates,
    sort_bam,
)
from exomeflow.cnv import run_cnv_calling
from exomeflow.fastp import run_fastp
from exomeflow.filtering import run_cohort_filtration, run_variant_filtration
from exomeflow.hpo_annotation import run_hpo_annotation, run_hpo_annotation_cohort
from exomeflow.joint_genotyping import run_joint_genotyping
from exomeflow.logger import close_sample_logger, get_pipeline_logger, get_sample_logger
from exomeflow.recalibration import run_base_recalibration
from exomeflow.reporting import run_multiqc
from exomeflow.somatic import run_mutect2, run_somatic_filtration
from exomeflow.steps import CohortStep, SampleStep
from exomeflow.utils import Checkpoint, PipelineStepError, detect_samples
from exomeflow.variant_calling import run_haplotype_caller

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")


def _cohort_active(cfg: "Config") -> bool:
    """Joint genotyping only replaces per-sample output for germline mode."""
    return cfg.joint_genotyping and cfg.mode == "germline"


def _sample_is_complete(sample: str, cfg: "Config", checkpoint: Checkpoint) -> bool:
    """
    A sample is complete when every currently-applicable SampleStep has its
    own per-step checkpoint marked done — derived dynamically from the step
    list rather than a separate coarse "COMPLETE" flag.

    Bug found via audit: the old coarse flag (still set for backward-
    compatible checkpoint directories, but no longer consulted) meant (a)
    upgrading to a version with new steps (e.g. v1 output dirs missing the
    new hpo/acmg steps) permanently skipped those samples instead of
    reprocessing just the new steps, and (b) a step that gracefully skipped
    this run (e.g. ACMG with InterVar not yet provisioned) got the sample
    marked done anyway, silently blocking any future retry once the
    prerequisite was fixed. Deriving completeness from the actual step list
    fixes both: only the steps genuinely missing their checkpoint re-run —
    everything already done is skipped internally by each step's own
    `if checkpoint.done(...): return` guard, so this doesn't mean starting
    a "complete" sample over from scratch.
    """
    return all(
        checkpoint.done(sample, step.name)
        for step in SAMPLE_STEPS
        if step.applies(cfg)
    )


SAMPLE_STEPS: list[SampleStep] = [
    SampleStep("fastp", run_fastp),
    SampleStep("bwa", run_bwa_mem),
    SampleStep("sort", sort_bam),
    SampleStep("flagstat", generate_flagstat),
    SampleStep("markdup", mark_duplicates),
    SampleStep("index", build_bam_index),
    SampleStep("bqsr", run_base_recalibration),
    SampleStep("haplotype", run_haplotype_caller, applies=lambda cfg: cfg.mode == "germline"),
    SampleStep("mutect2", run_mutect2, applies=lambda cfg: cfg.mode == "somatic"),
    SampleStep("cnv", run_cnv_calling, applies=lambda cfg: cfg.call_cnv),
    SampleStep(
        # Equivalent to "germline and not _cohort_active(cfg)" — simplified
        # since _cohort_active already implies mode == "germline", so the
        # mode check doesn't need restating here too (found via audit).
        "filter", run_variant_filtration,
        applies=lambda cfg: cfg.mode == "germline" and not cfg.joint_genotyping,
    ),
    SampleStep("somatic_filter", run_somatic_filtration, applies=lambda cfg: cfg.mode == "somatic"),
    SampleStep("annovar", run_annovar_annotation, applies=lambda cfg: not _cohort_active(cfg)),
    SampleStep("hpo", run_hpo_annotation, applies=lambda cfg: not _cohort_active(cfg)),
    SampleStep("acmg", run_intervar, applies=lambda cfg: not _cohort_active(cfg)),
]

COHORT_STEPS: list[CohortStep] = [
    CohortStep("joint_genotyping", run_joint_genotyping, applies=_cohort_active),
    CohortStep("cohort_filter", run_cohort_filtration, applies=_cohort_active),
    CohortStep("cohort_annovar", run_cohort_annotation, applies=_cohort_active),
    CohortStep("cohort_hpo", run_hpo_annotation_cohort, applies=_cohort_active),
    CohortStep("cohort_acmg", run_intervar_cohort, applies=_cohort_active),
    CohortStep("multiqc", run_multiqc),  # always attempted; best-effort, never fails the run
]


# ---------------------------------------------------------------------------
# Single-sample entry point
# ---------------------------------------------------------------------------

def process_sample(sample: str, cfg: "Config", timestamp: str) -> None:
    """
    Run every applicable SampleStep for a single *sample*.

    This function is designed to be safe to call inside a subprocess worker
    (ProcessPoolExecutor) — it re-initialises logging so each worker writes
    to its own log file.
    """
    # Re-init logging inside the worker process
    get_pipeline_logger(cfg.log_dir, timestamp)
    get_sample_logger(sample, cfg.log_dir, timestamp)

    sample_logger = logging.getLogger(f"exomeflow.{sample}")
    checkpoint = Checkpoint(cfg.checkpoint_dir, cfg.genome_build)

    try:
        sample_logger.info("=" * 50)
        sample_logger.info("Processing sample: %s (mode=%s)", sample, cfg.mode)
        sample_logger.info("=" * 50)

        for step in SAMPLE_STEPS:
            if not step.applies(cfg):
                continue
            try:
                step.run(sample, cfg, checkpoint)
            except PipelineStepError as exc:
                sample_logger.error("Failed at step '%s': %s", step.name, exc)
                raise

        sample_logger.log(25, "Sample %s completed successfully.", sample)
        sample_logger.info("=" * 10 + " OUTPUT FILES: %s " + "=" * 10, sample)
        sample_logger.info("  BAM (IGV):     %s", cfg.map_dir / f"{sample}_recalibrated.bam")
        if not _cohort_active(cfg):
            sample_logger.info("  PASS VCF:      %s", cfg.vcf_dir / f"{sample}_PASS.vcf")
            sample_logger.info(
                "  Annotated TXT: %s",
                cfg.vcf_dir / f"{sample}.annovar.{cfg.annovar_buildver}_multianno.txt",
            )
            sample_logger.info("  HPO/ACMG TXT:  %s", cfg.vcf_dir / f"{sample}.annovar.hpo.txt")
        sample_logger.info("=" * 50)
    finally:
        # A ProcessPoolExecutor worker is reused across many samples
        # whenever len(samples) > max_workers — without this, each new
        # sample's FileHandler (opened by get_sample_logger above) is never
        # released, leaking one file descriptor per sample the worker ever
        # processes over its lifetime. Found via audit.
        close_sample_logger(sample)


# ---------------------------------------------------------------------------
# Cohort orchestrator
# ---------------------------------------------------------------------------

def run_pipeline(cfg: "Config") -> int:
    """
    Run the full WES pipeline for all samples found in *cfg.input_dir*.

    Returns the number of failed samples (0 = success).
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    cfg.setup_directories()
    get_pipeline_logger(cfg.log_dir, timestamp)

    logger.info("=" * 50)
    logger.info("ExomeFlow — Automated Whole Exome Sequencing Analysis Python Package")
    logger.info("=" * 50)
    logger.info("Author  : Robin Kumar")
    logger.info("Affil.  : Dr Prabudh Goel Lab, All India Institute of Medical Sciences (AIIMS), New Delhi")
    logger.info("Email   : itsrobintomar@gmail.com")
    logger.info("Version : %s", __import__("exomeflow").__version__)
    logger.info("Mode    : %s | build=%s | joint-genotyping=%s | cnv=%s",
                cfg.mode, cfg.genome_build, cfg.joint_genotyping, cfg.call_cnv)
    logger.info("PyPI    : https://pypi.org/project/exomeflow")
    logger.info("=" * 50)
    logger.info("Parallel mode - up to %d samples simultaneously", cfg.max_workers)
    logger.info("Started at: %s", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    # Dependency/reference resolution already happened in cli.py's
    # _ensure_ready() -> setup_env.check_and_fix_dependencies() before
    # run_pipeline() was ever called — no second check here.

    # Sample discovery
    logger.info("Scanning %s for samples ...", cfg.input_dir)
    samples = detect_samples(cfg.input_dir)
    logger.info("Found %d sample(s): %s", len(samples), ", ".join(samples))
    logger.info("Will process %d sample(s) in parallel", cfg.max_workers)

    checkpoint = Checkpoint(cfg.checkpoint_dir, cfg.genome_build)

    # Skip already-completed samples
    pending = [s for s in samples if not _sample_is_complete(s, cfg, checkpoint)]
    skipped = len(samples) - len(pending)
    if skipped:
        logger.info("%d sample(s) already completed; skipping.", skipped)

    # --------------------------------------------------------------- dispatch
    failed: list[str] = []
    processed: list[str] = []

    # Graceful interrupt — registered unconditionally (not just when there's
    # per-sample work) since the cohort phase below (GenomicsDBImport,
    # ANNOVAR, InterVar, MultiQC) can also run for a long time even when
    # every sample was already complete on entry. Found via audit.
    def _handle_interrupt(sig, frame):  # type: ignore[no-untyped-def]
        logger.error("Pipeline interrupted by user.")
        sys.exit(1)

    signal.signal(signal.SIGINT, _handle_interrupt)
    signal.signal(signal.SIGTERM, _handle_interrupt)

    if pending:
        if cfg.max_workers == 1:
            # Sequential - run in the current process (simpler traceback)
            for sample in pending:
                try:
                    process_sample(sample, cfg, timestamp)
                    processed.append(sample)
                except Exception as exc:
                    logger.error("Sample %s failed: %s", sample, exc)
                    failed.append(sample)
        else:
            # Parallel - each sample runs in its own worker process
            with ProcessPoolExecutor(max_workers=cfg.max_workers) as pool:
                future_to_sample = {
                    pool.submit(process_sample, s, cfg, timestamp): s
                    for s in pending
                }
                for future in as_completed(future_to_sample):
                    sample = future_to_sample[future]
                    exc = future.exception()
                    if exc:
                        logger.error("Sample %s failed: %s", sample, exc)
                        failed.append(sample)
                    else:
                        logger.log(25, "Sample %s completed.", sample)
                        processed.append(sample)

    # ----------------------------------------------------------- cohort phase
    # Runs once after every sample's SampleSteps finish, over every sample
    # that didn't hard-fail this run (previously-completed + newly-processed,
    # excluding anything that just failed).
    #
    # Deliberately *not* `_sample_is_complete()` here (bug found via audit,
    # live smoke test): that check is strict - every applicable step,
    # including optional ones like ACMG/HPO, must be checkpointed. ACMG
    # gracefully skipping (e.g. InterVar not yet provisioned) doesn't raise
    # PipelineStepError and isn't a real failure, but it would make
    # _sample_is_complete() false for every sample, silently skipping the
    # *entire* cohort phase - including MultiQC - on an otherwise fully
    # successful run.
    cohort_samples = [s for s in samples if s not in failed]
    if cohort_samples:
        for cstep in COHORT_STEPS:
            if not cstep.applies(cfg):
                continue
            if checkpoint.done("__cohort__", cstep.name):
                logger.info("[cohort] %s already completed, skipping.", cstep.name)
                continue
            try:
                result = cstep.run(cohort_samples, cfg)
                # Steps that don't distinguish "succeeded" from "gracefully
                # skipped" via a return value implicitly return None here,
                # which is deliberately treated as success (their contract
                # is: raise PipelineStepError on real failure). Only a step
                # that explicitly returns False (e.g. MultiQC skipping
                # because the tool isn't installed) is left unmarked, so a
                # later retry once the blocker is fixed isn't silently
                # skipped forever.
                if result is not False:
                    checkpoint.mark("__cohort__", cstep.name)
            except PipelineStepError as exc:
                logger.error("[cohort] Failed at step '%s': %s", cstep.name, exc)
                failed.append(f"__cohort__:{cstep.name}")

    # --------------------------------------------------------------- summary
    logger.info("=" * 50)
    logger.info("Pipeline Summary")
    logger.info("=" * 50)
    logger.info("Total samples   : %d", len(samples))
    logger.info("Already done    : %d", skipped)
    logger.info("Successfully run: %d", len(processed))
    logger.info("Failed          : %d", len(failed))
    logger.info("Completed at    : %s", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("=" * 50)

    if failed:
        logger.error(
            "Pipeline completed with %d error(s). "
            "Check individual sample logs in: %s",
            len(failed), cfg.log_dir,
        )
        for s in failed:
            logger.error("  FAILED: %s", s)
    else:
        logger.log(25, "Pipeline completed successfully!")
        logger.info(
            "Output: %s",
            cfg.cohort_dir if _cohort_active(cfg) else cfg.vcf_dir,
        )

    return len(failed)
