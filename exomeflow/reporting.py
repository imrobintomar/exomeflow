"""
Cohort step — MultiQC rollup across fastp/flagstat/GATK metrics.

Always attempted (not gated behind a flag) but never fails the pipeline:
missing/failed multiqc only produces a warning, since this is a QC
convenience output, not a correctness-critical one.
"""

from __future__ import annotations

import logging
import shutil
import subprocess

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")

STEP = "multiqc"


def run_multiqc(samples: list[str], cfg: "Config") -> None:
    """Aggregate all sample QC/log outputs into one MultiQC HTML report."""
    if not shutil.which("multiqc"):
        logger.warning(
            "[cohort] multiqc not found on PATH — skipping rollup report "
            "(non-fatal; run `pip install multiqc` to enable it)."
        )
        return

    report_dir = cfg.output_dir / "multiqc"
    report_dir.mkdir(parents=True, exist_ok=True)

    logger.info("[cohort] Running MultiQC over %s ...", cfg.output_dir)

    result = subprocess.run(
        [
            "multiqc", str(cfg.output_dir),
            "-o", str(report_dir),
            "-n", "exomeflow_report",
            "--force",
        ],
        env=cfg.env(),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        logger.warning(
            "[cohort] multiqc exited non-zero (%d) — skipping rollup report:\n%s",
            result.returncode, result.stderr[-500:],
        )
        return

    logger.log(25, "[cohort] MultiQC report: %s", report_dir / "exomeflow_report.html")
