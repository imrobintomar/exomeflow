"""
Shared utilities: subprocess wrapper, sample detection, dependency checks,
checkpointing, and the full requirements check run as pipeline step 0.
"""

from __future__ import annotations

import importlib
import logging
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")


# ---------------------------------------------------------------------------
# Subprocess wrapper
# ---------------------------------------------------------------------------

class PipelineStepError(RuntimeError):
    """Raised when an external tool returns a non-zero exit code."""


def run_cmd(
    cmd: list[str],
    *,
    env: dict[str, str] | None = None,
    cwd: Path | None = None,
    step_name: str = "",
    sample: str = "",
) -> subprocess.CompletedProcess:  # type: ignore[type-arg]
    """
    Execute *cmd* via subprocess.run(), streaming stdout/stderr to the
    caller's inherited file descriptors.

    Raises PipelineStepError on non-zero exit so callers can catch a
    specific exception rather than inspecting return codes.
    """
    label = f"[{sample}] {step_name}" if sample else step_name
    logger.debug("%s: %s", label, " ".join(str(c) for c in cmd))

    result = subprocess.run(
        cmd,
        env=env,
        cwd=str(cwd) if cwd else None,
        check=False,          # we check manually below
    )

    if result.returncode != 0:
        raise PipelineStepError(
            f"{label} failed with exit code {result.returncode}. "
            f"Command: {' '.join(str(c) for c in cmd)}"
        )
    return result


# ---------------------------------------------------------------------------
# Sample detection
# ---------------------------------------------------------------------------

def detect_samples(input_dir: Path) -> list[str]:
    """
    Scan *input_dir* for paired FASTQ files and return sorted sample IDs.

    Expected naming convention::

        <sample_id>_1.fastq.gz
        <sample_id>_2.fastq.gz

    Any base name matching ``*_1.fastq.gz`` is considered a sample; its
    partner ``*_2.fastq.gz`` is expected to exist (validated later per-step).
    """
    pattern = re.compile(r"^(.+)_1\.fastq\.gz$")
    samples: set[str] = set()

    for f in input_dir.glob("*.fastq.gz"):
        m = pattern.match(f.name)
        if m:
            samples.add(m.group(1))

    if not samples:
        raise FileNotFoundError(
            f"No paired FASTQ files (matching *_1.fastq.gz) found in {input_dir}"
        )

    return sorted(samples)


# ---------------------------------------------------------------------------
# Dependency validation
# ---------------------------------------------------------------------------

def check_dependencies() -> None:
    """Verify that all required external tools are on PATH."""
    required = ["bwa", "samtools", "gatk", "fastp"]
    missing = [cmd for cmd in required if not _which(cmd)]

    if missing:
        raise EnvironmentError(
            "The following required tools were not found on PATH: "
            + ", ".join(missing)
        )
    logger.info("All dependencies found: %s", ", ".join(required))


def _which(cmd: str) -> bool:
    import shutil
    return shutil.which(cmd) is not None


def check_reference_files(cfg: "Config") -> None:
    """Verify that all reference files and directories exist."""
    to_check: list[Path] = [
        cfg.reference,
        cfg.dbsnp,
        cfg.mills,
        cfg.known_indels,
        cfg.annovar_db,
    ]
    missing = [p for p in to_check if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Reference files not found:\n"
            + "\n".join(f"  {p}" for p in missing)
        )
    logger.info("All reference files validated.")


# ---------------------------------------------------------------------------
# Full requirements check (runs as pipeline step 0)
# ---------------------------------------------------------------------------

# Minimum versions for system tools
_SYSTEM_TOOLS: list[tuple[str, list[str], str]] = [
    # (name, version_cmd, min_version)
    ("bwa",      ["bwa"],                "0.7.17"),
    ("samtools", ["samtools", "--version"], "1.13"),
    ("gatk",     ["gatk", "--version"],  "4.6.0"),
    ("fastp",    ["fastp", "--version"], "0.20.1"),
    ("perl",     ["perl", "--version"],  "5.26"),
]

# Minimum versions for Python packages
_PYTHON_PACKAGES: list[tuple[str, str]] = [
    # (import_name, min_version)
    ("typer",      "0.12.0"),
    ("rich",       "13.0.0"),
    ("pandas",     "2.0.0"),
    ("matplotlib", "3.7.0"),
]


def _parse_version(text: str) -> tuple[int, ...]:
    m = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", text)
    if not m:
        return (0,)
    return tuple(int(g) for g in m.groups() if g is not None)


def _version_ok(found: tuple[int, ...], minimum: str) -> bool:
    return found >= _parse_version(minimum)


def check_all_requirements(cfg: "Config") -> None:
    """
    Run a full pre-flight requirements check before the pipeline starts.

    Checks:
      - Python version >= 3.9
      - All required Python packages (with version)
      - All system tools on PATH (with version)
      - ANNOVAR installation directory
      - Reference files

    Raises EnvironmentError listing every failure found so the user can
    fix them all at once rather than one at a time.
    """
    failures: list[str] = []
    logger.info("=" * 50)
    logger.info("Step 0 — Requirements Check")
    logger.info("=" * 50)

    # ── Python version ────────────────────────────────────────────────────
    major, minor, micro = sys.version_info[:3]
    py_ver = f"{major}.{minor}.{micro}"
    if (major, minor) < (3, 9):
        failures.append(
            f"Python {py_ver} is below minimum 3.9  "
            f"→  conda install python>=3.9"
        )
        logger.error("Python %s < 3.9", py_ver)
    else:
        logger.info("  ✔  Python %s", py_ver)

    # ── Python packages ───────────────────────────────────────────────────
    for import_name, min_ver in _PYTHON_PACKAGES:
        try:
            mod = importlib.import_module(import_name)
            version_str = getattr(mod, "__version__", None)
            if version_str is None:
                from importlib.metadata import version as _meta_ver
                version_str = _meta_ver(import_name)
            found = _parse_version(version_str or "0")
            if not _version_ok(found, min_ver):
                msg = (
                    f"{import_name} {version_str} < {min_ver}  "
                    f"→  pip install --upgrade {import_name}>={min_ver}"
                )
                failures.append(msg)
                logger.error("  ✘  %s", msg)
            else:
                logger.info("  ✔  %s %s", import_name, version_str)
        except ImportError:
            msg = (
                f"{import_name} is not installed  "
                f"→  pip install {import_name}>={min_ver}"
            )
            failures.append(msg)
            logger.error("  ✘  %s", msg)

    # ── System tools ──────────────────────────────────────────────────────
    _INSTALL_HINTS = {
        "bwa":      "conda install -c bioconda bwa",
        "samtools": "conda install -c bioconda samtools",
        "gatk":     "https://github.com/broadinstitute/gatk/releases  then add to PATH",
        "fastp":    "conda install -c bioconda fastp",
        "perl":     "conda install perl",
    }

    for name, version_cmd, min_ver in _SYSTEM_TOOLS:
        if not shutil.which(name):
            msg = (
                f"{name} not found on PATH  "
                f"→  {_INSTALL_HINTS.get(name, 'install ' + name)}"
            )
            failures.append(msg)
            logger.error("  ✘  %s", msg)
            continue

        try:
            result = subprocess.run(
                version_cmd,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=10,
            )
            output = result.stdout + result.stderr
            found = _parse_version(output)
            ver_str = ".".join(str(v) for v in found)
            if not _version_ok(found, min_ver):
                msg = (
                    f"{name} {ver_str} < {min_ver}  "
                    f"→  {_INSTALL_HINTS.get(name, 'update ' + name)}"
                )
                failures.append(msg)
                logger.error("  ✘  %s", msg)
            else:
                logger.info("  ✔  %s %s", name, ver_str)
        except Exception as exc:
            msg = f"{name} found but version check failed: {exc}"
            failures.append(msg)
            logger.error("  ✘  %s", msg)

    # ── ANNOVAR ───────────────────────────────────────────────────────────
    annovar_pl = Path(cfg.annovar_bin) / "table_annovar.pl"
    if not annovar_pl.exists():
        msg = (
            f"ANNOVAR not found at {annovar_pl}  "
            f"→  check --annovar-bin path or download from "
            f"https://annovar.openbioinformatics.org"
        )
        failures.append(msg)
        logger.error("  ✘  %s", msg)
    else:
        logger.info("  ✔  ANNOVAR found at %s", annovar_pl)

    # ── Reference files ───────────────────────────────────────────────────
    ref_files = {
        "Reference genome (--reference)": cfg.reference,
        "dbSNP VCF (--dbsnp)":            cfg.dbsnp,
        "Mills indels (--mills)":          cfg.mills,
        "Known indels (--known-indels)":   cfg.known_indels,
        "ANNOVAR humandb (--annovar-db)":  cfg.annovar_db,
    }
    for label, path in ref_files.items():
        if not Path(path).exists():
            msg = f"{label} not found: {path}"
            failures.append(msg)
            logger.error("  ✘  %s", msg)
        else:
            logger.info("  ✔  %s: %s", label, path)

    # ── Exome intervals (optional — warn only) ───────────────────────────
    if cfg.intervals and not Path(cfg.intervals).exists():
        logger.warning(
            "  ⚠  Exome intervals BED not found: %s — "
            "pipeline will run in whole-genome mode (slower).",
            cfg.intervals,
        )

    # ── Result ────────────────────────────────────────────────────────────
    if failures:
        summary = "\n".join(f"  • {f}" for f in failures)
        raise EnvironmentError(
            f"Requirements check failed — fix the following issues and retry:\n"
            f"{summary}"
        )

    logger.log(25, "Step 0 — All requirements satisfied. Starting pipeline.")


# ---------------------------------------------------------------------------
# Checkpointing
# ---------------------------------------------------------------------------

class Checkpoint:
    """
    Lightweight file-based checkpoint system.

    A step is considered complete when the file
    ``<checkpoint_dir>/<sample>.<step>.done`` exists.
    """

    def __init__(self, checkpoint_dir: Path) -> None:
        self._dir = checkpoint_dir
        self._dir.mkdir(parents=True, exist_ok=True)

    def _path(self, sample: str, step: str) -> Path:
        return self._dir / f"{sample}.{step}.done"

    def mark(self, sample: str, step: str) -> None:
        """Create the checkpoint file for *sample* / *step*."""
        self._path(sample, step).touch()

    def done(self, sample: str, step: str) -> bool:
        """Return True if the checkpoint file exists."""
        return self._path(sample, step).exists()

    def mark_sample_complete(self, sample: str) -> None:
        self.mark(sample, "COMPLETE")

    def is_sample_complete(self, sample: str) -> bool:
        return self.done(sample, "COMPLETE")


# ---------------------------------------------------------------------------
# Variant count helper (used in filtering summary)
# ---------------------------------------------------------------------------

def count_variants(vcf: Path) -> int:
    """Return the number of non-header lines in *vcf* (0 if file missing)."""
    if not vcf.exists():
        return 0
    count = 0
    with open(vcf, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("#"):
                count += 1
    return count
