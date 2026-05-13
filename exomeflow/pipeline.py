"""
Orchestration layer — processes a single sample end-to-end and supports
parallel cohort processing via concurrent.futures.ProcessPoolExecutor.

Design mirrors the Bash ``process_sample`` / ``main`` functions:
  - Each sample runs all 12 steps in order.
  - A failed step raises PipelineStepError; the exception is caught at the
    cohort level so other samples continue processing.
  - Checkpointing lets the pipeline resume where it left off.
"""

from __future__ import annotations

import logging
import signal
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from typing import TYPE_CHECKING

from exomeflow.annotation import run_annovar_annotation
from exomeflow.bam_processing import (
    build_bam_index,
    generate_flagstat,
    mark_duplicates,
    sort_bam,
)
from exomeflow.alignment import run_bwa_mem
from exomeflow.fastp import run_fastp
from exomeflow.filtering import run_variant_filtration
from exomeflow.logger import get_pipeline_logger, get_sample_logger
from exomeflow.recalibration import run_base_recalibration
from exomeflow.utils import (
    Checkpoint,
    PipelineStepError,
    check_all_requirements,
    detect_samples,
)
from exomeflow.variant_calling import run_haplotype_caller

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")


# ---------------------------------------------------------------------------
# Single-sample entry point
# ---------------------------------------------------------------------------

def process_sample(sample: str, cfg: "Config", timestamp: str) -> None:
    """
    Run the complete WES pipeline for a single *sample*.

    This function is designed to be safe to call inside a subprocess worker
    (ProcessPoolExecutor) — it re-initialises logging so each worker writes
    to its own log file.
    """
    # Re-init logging inside the worker process
    get_pipeline_logger(cfg.log_dir, timestamp)
    get_sample_logger(sample, cfg.log_dir, timestamp)

    sample_logger = logging.getLogger(f"exomeflow.{sample}")
    checkpoint = Checkpoint(cfg.checkpoint_dir)

    sample_logger.info("=" * 50)
    sample_logger.info("Processing sample: %s", sample)
    sample_logger.info("=" * 50)

    steps = [
        ("fastp", lambda: run_fastp(sample, cfg, checkpoint)),
        ("bwa", lambda: run_bwa_mem(sample, cfg, checkpoint)),
        ("sort", lambda: sort_bam(sample, cfg, checkpoint)),
        ("flagstat", lambda: generate_flagstat(sample, cfg, checkpoint)),
        ("markdup", lambda: mark_duplicates(sample, cfg, checkpoint)),
        ("index", lambda: build_bam_index(sample, cfg, checkpoint)),
        ("bqsr", lambda: run_base_recalibration(sample, cfg, checkpoint)),
        ("haplotype", lambda: run_haplotype_caller(sample, cfg, checkpoint)),
        ("filter", lambda: run_variant_filtration(sample, cfg, checkpoint)),
        ("annovar", lambda: run_annovar_annotation(sample, cfg, checkpoint)),
    ]

    for step_name, step_fn in steps:
        try:
            step_fn()
        except PipelineStepError as exc:
            sample_logger.error("Failed at step '%s': %s", step_name, exc)
            raise

    checkpoint.mark_sample_complete(sample)

    sample_logger.log(25, "Sample %s completed successfully.", sample)
    sample_logger.info("=" * 10 + " OUTPUT FILES: %s " + "=" * 10, sample)
    sample_logger.info("  BAM (IGV):     %s", cfg.map_dir / f"{sample}_recalibrated.bam")
    sample_logger.info("  BAM index:     %s", cfg.map_dir / f"{sample}_recalibrated.bam.bai")
    sample_logger.info("  Raw VCF:       %s", cfg.vcf_dir / f"{sample}.vcf")
    sample_logger.info("  PASS VCF:      %s", cfg.vcf_dir / f"{sample}_PASS.vcf")
    sample_logger.info("  Annotated VCF: %s", cfg.vcf_dir / f"{sample}.annovar.hg38_multianno.vcf")
    sample_logger.info("  Annotated TXT: %s", cfg.vcf_dir / f"{sample}.annovar.hg38_multianno.txt")
    sample_logger.info("=" * 50)


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
    logger.info("ExomeFlow — Whole Exome Sequencing Analysis Pipeline")
    logger.info("=" * 50)
    logger.info("Author  : Robin Kumar")
    logger.info("Affil.  : All India Institute of Medical Sciences (AIIMS), New Delhi")
    logger.info("Email   : itsrobintomar@gmail.com")
    logger.info("Version : %s", __import__("exomeflow").__version__)
    logger.info("Citation: Tomar R. ExomeFlow: A Python-based production-quality")
    logger.info("          Whole Exome Sequencing analysis pipeline.")
    logger.info("          GitHub: https://github.com/imrobintomar/exomeflow")
    logger.info("          PyPI:   https://pypi.org/project/exomeflow")
    logger.info("=" * 50)
    logger.info("Parallel mode — up to %d samples simultaneously", cfg.max_workers)
    logger.info("Started at: %s", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    # Step 0 — full requirements check (tools, packages, reference files)
    check_all_requirements(cfg)

    # Sample discovery
    logger.info("Scanning %s for samples ...", cfg.input_dir)
    samples = detect_samples(cfg.input_dir)
    logger.info("Found %d sample(s): %s", len(samples), ", ".join(samples))
    logger.info("Will process %d sample(s) in parallel", cfg.max_workers)

    checkpoint = Checkpoint(cfg.checkpoint_dir)

    # Skip already-completed samples
    pending = [s for s in samples if not checkpoint.is_sample_complete(s)]
    skipped = len(samples) - len(pending)
    if skipped:
        logger.info("%d sample(s) already completed; skipping.", skipped)

    if not pending:
        logger.log(25, "All samples already completed.")
        return 0

    # Graceful interrupt
    def _handle_interrupt(sig, frame):  # type: ignore[no-untyped-def]
        logger.error("Pipeline interrupted by user.")
        sys.exit(1)

    signal.signal(signal.SIGINT,  _handle_interrupt)
    signal.signal(signal.SIGTERM, _handle_interrupt)

    # --------------------------------------------------------------- dispatch
    failed: list[str] = []
    processed: list[str] = []

    if cfg.max_workers == 1:
        # Sequential — run in the current process (simpler traceback)
        for sample in pending:
            try:
                process_sample(sample, cfg, timestamp)
                processed.append(sample)
            except Exception as exc:
                logger.error("Sample %s failed: %s", sample, exc)
                failed.append(sample)
    else:
        # Parallel — each sample runs in its own worker process
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
        logger.info("Annotated VCF files are in: %s", cfg.vcf_dir)

    return len(failed)
