"""
Shared utilities: subprocess wrapper, sample detection, checkpointing, and
version-parsing helpers used by setup_env's dependency checker.
"""

from __future__ import annotations

import gzip
import logging
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger("exomeflow")


# ---------------------------------------------------------------------------
# System resource detection — used to size --threads/--java-opts sensibly
# instead of a one-size-fits-all static default that can undersubscribe a
# big box or, worse, request more RAM than a small one actually has.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SystemResources:
    cpu_count: int
    total_ram_gb: float
    available_ram_gb: float
    free_disk_gb: float


def detect_system_resources(disk_path: Path = Path(".")) -> SystemResources:
    """Read CPU/RAM/disk from the OS. Linux-only (matches this package's scope)."""
    cpu_count = os.cpu_count() or 1

    total_ram_gb = available_ram_gb = 0.0
    try:
        meminfo = {}
        with open("/proc/meminfo", encoding="utf-8") as fh:
            for line in fh:
                key, _, rest = line.partition(":")
                if key in ("MemTotal", "MemAvailable"):
                    meminfo[key] = int(rest.strip().split()[0])  # kB
        total_ram_gb = meminfo.get("MemTotal", 0) / (1024 * 1024)
        available_ram_gb = meminfo.get("MemAvailable", 0) / (1024 * 1024)
    except FileNotFoundError:
        pass  # non-Linux — leave at 0.0, callers fall back to conservative defaults

    # disk_usage() requires an existing path — --output is documented as
    # "created if absent" and usually doesn't exist yet at this point, so
    # walk up to the nearest existing ancestor instead of crashing outright.
    probe = Path(disk_path).absolute()
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    free_disk_gb = shutil.disk_usage(probe).free / (1024 ** 3)

    return SystemResources(cpu_count, total_ram_gb, available_ram_gb, free_disk_gb)


def recommend_threads(resources: SystemResources) -> int:
    """Leave 2 cores free for the OS/other processes; cap at 24 (diminishing
    returns for BWA/HaplotypeCaller beyond this in practice)."""
    return max(1, min(resources.cpu_count - 2, 24))


def recommend_java_opts(resources: SystemResources) -> str:
    """
    Use 60% of *available* (not total) RAM, leaving headroom for BWA/samtools
    running alongside the JVM. Floors at 4g, caps at 80g.
    """
    if resources.available_ram_gb <= 0:
        return "-Xmx8g"  # unknown RAM (non-Linux) — conservative fallback
    gb = max(4, min(int(resources.available_ram_gb * 0.6), 80))
    return f"-Xmx{gb}g"


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
# Version helpers (shared by setup_env's dependency checker)
# ---------------------------------------------------------------------------

def _parse_version(text: str) -> tuple[int, ...]:
    m = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", text)
    if not m:
        return (0,)
    return tuple(int(g) for g in m.groups() if g is not None)


def _version_ok(found: tuple[int, ...], minimum: str) -> bool:
    return found >= _parse_version(minimum)


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


# ---------------------------------------------------------------------------
# Variant count helper (used in filtering summary)
# ---------------------------------------------------------------------------

def count_variants(vcf: Path) -> int:
    """
    Return the number of non-header lines in *vcf* (0 if file missing).

    Handles both plain-text and bgzip/gzip VCFs (e.g. GATK auto-bgzips any
    -O output whose filename ends in .gz, such as the cohort joint-genotyped
    VCF) — a plain text open() on a gzipped file raises UnicodeDecodeError.
    """
    if not vcf.exists():
        return 0
    opener = gzip.open if vcf.suffix == ".gz" else open
    count = 0
    with opener(vcf, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("#"):
                count += 1
    return count
