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
# Full requirements check with auto-install (runs as pipeline step 0)
# ---------------------------------------------------------------------------

# Minimum versions for system tools
# (name, version_cmd, min_version, conda_pkg)
_SYSTEM_TOOLS: list[tuple[str, list[str], str, str]] = [
    ("bwa",      ["bwa"],                  "0.7.17", "bioconda::bwa"),
    ("samtools", ["samtools", "--version"], "1.13",   "bioconda::samtools"),
    ("fastp",    ["fastp", "--version"],   "0.20.1", "bioconda::fastp"),
    ("perl",     ["perl", "--version"],    "5.26",   "conda-forge::perl"),
    # gatk is not conda-installable here — requires manual PATH setup
    ("gatk",     ["gatk", "--version"],    "4.6.0",  ""),
]

# Python packages (import_name, pypi_name, min_version)
_PYTHON_PACKAGES: list[tuple[str, str, str]] = [
    ("typer",      "typer",      "0.12.0"),
    ("rich",       "rich",       "13.0.0"),
    ("pandas",     "pandas",     "2.0.0"),
    ("matplotlib", "matplotlib", "3.7.0"),
]

# ANNOVAR databases to check inside humandb/
# (filename_pattern, db_name, description)
_ANNOVAR_DATABASES: list[tuple[str, str, str]] = [
    ("hg38_refGene.txt",             "refGene",            "Gene annotation"),
    ("hg38_avsnp150.txt",            "avsnp150",           "dbSNP 150"),
    ("hg38_clinvar_20240416.txt",    "clinvar_20240416",   "ClinVar"),
    ("hg38_gnomad41_exome.txt",      "gnomad41_exome",     "gnomAD v4.1 exome"),
    ("hg38_gnomad41_genome.txt",     "gnomad41_genome",    "gnomAD v4.1 genome"),
    ("hg38_dbnsfp47a.txt",           "dbnsfp47a",          "dbNSFP 4.7a"),
    ("hg38_cosmic84_coding.txt",     "cosmic84_coding",    "COSMIC v84"),
    ("hg38_exac03.txt",              "exac03",             "ExAC 0.3"),
]


def _parse_version(text: str) -> tuple[int, ...]:
    m = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", text)
    if not m:
        return (0,)
    return tuple(int(g) for g in m.groups() if g is not None)


def _version_ok(found: tuple[int, ...], minimum: str) -> bool:
    return found >= _parse_version(minimum)


def _ask(question: str) -> bool:
    """Print *question* and return True if the user answers y/yes."""
    try:
        answer = input(f"\n  {question} [y/N]: ").strip().lower()
        return answer in ("y", "yes")
    except (EOFError, KeyboardInterrupt):
        return False


def _pip_install(package: str, min_ver: str) -> bool:
    """Install *package* via pip. Returns True on success."""
    spec = f"{package}>={min_ver}"
    logger.info("  → Installing %s via pip ...", spec)
    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", spec],
        capture_output=True, text=True,
    )
    if result.returncode == 0:
        logger.log(25, "  ✔  %s installed successfully.", package)
        return True
    logger.error("  ✘  pip install failed:\n%s", result.stderr[-500:])
    return False


def _conda_install(conda_pkg: str) -> bool:
    """Install *conda_pkg* via conda. Returns True on success."""
    conda = shutil.which("conda") or shutil.which("mamba")
    if not conda:
        logger.warning("  ⚠  conda/mamba not found — cannot auto-install.")
        return False
    logger.info("  → Installing %s via conda ...", conda_pkg)
    result = subprocess.run(
        [conda, "install", "-y", "-c", conda_pkg.split("::")[0],
         conda_pkg.split("::")[-1]],
        capture_output=True, text=True,
    )
    if result.returncode == 0:
        logger.log(25, "  ✔  %s installed successfully.", conda_pkg)
        return True
    logger.error("  ✘  conda install failed:\n%s", result.stderr[-500:])
    return False


def _download_annovar_db(db_name: str, annovar_bin: Path, annovar_db: Path) -> bool:
    """
    Download a single ANNOVAR database using annotate_variation.pl.
    Returns True on success.
    """
    annotate_pl = annovar_bin / "annotate_variation.pl"
    if not annotate_pl.exists():
        logger.error("  ✘  annotate_variation.pl not found at %s", annotate_pl)
        return False

    logger.info("  → Downloading ANNOVAR database: %s ...", db_name)
    result = subprocess.run(
        [
            "perl", str(annotate_pl),
            "-buildver", "hg38",
            "-downdb",
            "-webfrom", "annovar",
            db_name,
            str(annovar_db),
        ],
        text=True,
    )
    if result.returncode == 0:
        logger.log(25, "  ✔  %s downloaded successfully.", db_name)
        return True
    logger.error("  ✘  Download of %s failed.", db_name)
    return False


def check_all_requirements(cfg: "Config") -> None:
    """
    Run a full pre-flight requirements check before the pipeline starts.

    - Missing Python packages  → auto-installs via pip (no prompt needed)
    - Missing system tools     → asks permission, then installs via conda
    - Missing ANNOVAR databases → asks permission per database, then downloads
    - Missing reference files  → reports with download instructions (cannot
                                  auto-download due to licensing)

    Raises EnvironmentError if any hard requirement remains unresolved after
    auto-install attempts.
    """
    failures: list[str] = []

    logger.info("=" * 50)
    logger.info("Step 0 — Requirements Check")
    logger.info("=" * 50)

    # ── Python version ────────────────────────────────────────────────────
    major, minor, micro = sys.version_info[:3]
    py_ver = f"{major}.{minor}.{micro}"
    if (major, minor) < (3, 9):
        failures.append(f"Python {py_ver} < 3.9  →  conda install python>=3.9")
        logger.error("  ✘  Python %s < 3.9", py_ver)
    else:
        logger.info("  ✔  Python %s", py_ver)

    # ── Python packages (auto-install via pip) ────────────────────────────
    for import_name, pypi_name, min_ver in _PYTHON_PACKAGES:
        try:
            mod = importlib.import_module(import_name)
            version_str = getattr(mod, "__version__", None)
            if version_str is None:
                from importlib.metadata import version as _meta_ver
                version_str = _meta_ver(pypi_name)
            found = _parse_version(version_str or "0")
            if not _version_ok(found, min_ver):
                logger.warning(
                    "  ⚠  %s %s < %s — upgrading automatically ...",
                    pypi_name, version_str, min_ver,
                )
                if not _pip_install(pypi_name, min_ver):
                    failures.append(
                        f"{pypi_name} upgrade failed  "
                        f"→  pip install --upgrade {pypi_name}>={min_ver}"
                    )
            else:
                logger.info("  ✔  %s %s", pypi_name, version_str)
        except ImportError:
            logger.warning("  ⚠  %s not found — installing automatically ...", pypi_name)
            if not _pip_install(pypi_name, min_ver):
                failures.append(
                    f"{pypi_name} not installed  "
                    f"→  pip install {pypi_name}>={min_ver}"
                )

    # ── System tools (ask → conda install) ───────────────────────────────
    for name, version_cmd, min_ver, conda_pkg in _SYSTEM_TOOLS:
        if not shutil.which(name):
            if not conda_pkg:
                # GATK — cannot auto-install, needs manual PATH setup
                failures.append(
                    f"{name} not found on PATH  "
                    f"→  download from https://github.com/broadinstitute/gatk/releases"
                    f" and add to PATH"
                )
                logger.error(
                    "  ✘  %s not found on PATH — must be installed manually.", name
                )
                continue

            logger.warning("  ⚠  %s not found on PATH.", name)
            if _ask(f"Install {name} now via conda? ({conda_pkg})"):
                if _conda_install(conda_pkg):
                    # Reload PATH in current process is not possible;
                    # warn user to restart shell if needed
                    logger.warning(
                        "  ⚠  %s installed. If the pipeline fails, open a new "
                        "terminal so PATH refreshes, then re-run.", name,
                    )
                else:
                    failures.append(
                        f"{name} install failed  "
                        f"→  conda install -c {conda_pkg.split('::')[0]} "
                        f"{conda_pkg.split('::')[-1]}"
                    )
            else:
                failures.append(
                    f"{name} not installed  "
                    f"→  conda install -c {conda_pkg.split('::')[0]} "
                    f"{conda_pkg.split('::')[-1]}"
                )
            continue

        # Tool is on PATH — check version
        try:
            result = subprocess.run(
                version_cmd,
                capture_output=True, text=True,
                encoding="utf-8", errors="replace", timeout=10,
            )
            output = result.stdout + result.stderr
            found  = _parse_version(output)
            ver_str = ".".join(str(v) for v in found)
            if not _version_ok(found, min_ver):
                logger.warning(
                    "  ⚠  %s %s < %s (minimum).", name, ver_str, min_ver
                )
                if conda_pkg and _ask(
                    f"Update {name} to >= {min_ver} via conda? ({conda_pkg})"
                ):
                    if not _conda_install(conda_pkg):
                        failures.append(f"{name} update failed")
                else:
                    failures.append(f"{name} {ver_str} < {min_ver}")
            else:
                logger.info("  ✔  %s %s", name, ver_str)
        except Exception as exc:
            failures.append(f"{name} version check failed: {exc}")
            logger.error("  ✘  %s version check failed: %s", name, exc)

    # ── ANNOVAR binary ────────────────────────────────────────────────────
    annovar_pl = Path(cfg.annovar_bin) / "table_annovar.pl"
    if not annovar_pl.exists():
        failures.append(
            f"ANNOVAR not found at {annovar_pl}  "
            f"→  download from https://annovar.openbioinformatics.org "
            f"and set --annovar-bin"
        )
        logger.error("  ✘  ANNOVAR not found at %s", annovar_pl)
    else:
        logger.info("  ✔  ANNOVAR found at %s", annovar_pl)

        # ── ANNOVAR databases (ask per missing db) ────────────────────────
        annovar_db = Path(cfg.annovar_db)
        missing_dbs: list[tuple[str, str]] = []

        for filename, db_name, description in _ANNOVAR_DATABASES:
            db_file = annovar_db / filename
            if not db_file.exists():
                missing_dbs.append((db_name, description))
                logger.warning("  ⚠  ANNOVAR db missing: %s (%s)", db_name, description)
            else:
                logger.info("  ✔  ANNOVAR db: %s", db_name)

        if missing_dbs:
            logger.warning(
                "\n  %d ANNOVAR database(s) are missing:", len(missing_dbs)
            )
            for db_name, desc in missing_dbs:
                logger.warning("    • %s — %s", db_name, desc)

            if _ask(
                f"Download all {len(missing_dbs)} missing ANNOVAR database(s) now?\n"
                f"  (This may take several minutes and requires internet access)"
            ):
                for db_name, desc in missing_dbs:
                    logger.info("  Downloading: %s (%s) ...", db_name, desc)
                    if not _download_annovar_db(
                        db_name, Path(cfg.annovar_bin), annovar_db
                    ):
                        failures.append(
                            f"ANNOVAR db download failed: {db_name}  "
                            f"→  run manually: perl annotate_variation.pl "
                            f"-buildver hg38 -downdb -webfrom annovar "
                            f"{db_name} {annovar_db}"
                        )
            else:
                for db_name, desc in missing_dbs:
                    failures.append(
                        f"ANNOVAR db missing: {db_name} ({desc})  "
                        f"→  perl annotate_variation.pl -buildver hg38 "
                        f"-downdb -webfrom annovar {db_name} {annovar_db}"
                    )

    # ── Reference files (cannot auto-download — licensing) ───────────────
    ref_files = {
        "Reference genome (--reference)": (cfg.reference,
            "Download: gsutil cp gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta ."),
        "dbSNP VCF (--dbsnp)": (cfg.dbsnp,
            "Download: gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz ."),
        "Mills indels (--mills)": (cfg.mills,
            "Download: gsutil cp gs://gatk-best-practices/somatic-hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ."),
        "Known indels (--known-indels)": (cfg.known_indels,
            "Download: gsutil cp gs://gatk-best-practices/somatic-hg38/Homo_sapiens_assembly38.known_indels.vcf.gz ."),
        "ANNOVAR humandb (--annovar-db)": (cfg.annovar_db, ""),
    }
    for label, (path, hint) in ref_files.items():
        if not Path(path).exists():
            msg = f"{label} not found: {path}"
            if hint:
                msg += f"\n      {hint}"
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

    # ── Final result ──────────────────────────────────────────────────────
    if failures:
        summary = "\n".join(f"  • {f}" for f in failures)
        raise EnvironmentError(
            f"Requirements check failed — unresolved issues:\n{summary}"
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
