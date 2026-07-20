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


def _cgroup_memory_limit_gb() -> float | None:
    """
    Return the container's cgroup memory cap in GB, or None if unset/not in
    a container. `/proc/meminfo` reports the *host's* RAM even inside a
    memory-limited Docker/Kubernetes container, so relying on it alone can
    recommend a JVM heap far beyond what the container will actually be
    allowed to use, leading to an OOM-killed run hours in. Found via audit.
    """
    for path in (Path("/sys/fs/cgroup/memory.max"),                    # cgroup v2
                 Path("/sys/fs/cgroup/memory/memory.limit_in_bytes")):  # cgroup v1
        try:
            raw = path.read_text().strip()
        except (FileNotFoundError, PermissionError, OSError):
            continue
        if raw == "max":
            continue  # v2's explicit "no limit"
        try:
            limit_bytes = int(raw)
        except ValueError:
            continue
        if limit_bytes > 2**62:
            continue  # v1's "no limit" sentinel (close to INT64_MAX)
        return limit_bytes / (1024 ** 3)
    return None


def detect_system_resources(disk_path: Path = Path(".")) -> SystemResources:
    """Read CPU/RAM/disk from the OS. Linux-only (matches this package's scope)."""
    try:
        # cgroup/affinity-aware — os.cpu_count() reports the host's full
        # core count even inside a `docker run --cpus=N`-limited container.
        # Found via audit.
        cpu_count = len(os.sched_getaffinity(0)) or 1
    except AttributeError:
        cpu_count = os.cpu_count() or 1  # non-Linux fallback

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

    cgroup_limit_gb = _cgroup_memory_limit_gb()
    if cgroup_limit_gb is not None:
        total_ram_gb = min(total_ram_gb, cgroup_limit_gb) if total_ram_gb else cgroup_limit_gb
        available_ram_gb = min(available_ram_gb, cgroup_limit_gb) if available_ram_gb else cgroup_limit_gb

    # disk_usage() requires an existing path — --output is documented as
    # "created if absent" and usually doesn't exist yet at this point, so
    # walk up to the nearest existing ancestor instead of crashing outright.
    probe = Path(disk_path).absolute()
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    free_disk_gb = shutil.disk_usage(probe).free / (1024 ** 3)

    return SystemResources(cpu_count, total_ram_gb, available_ram_gb, free_disk_gb)


def recommend_threads(resources: SystemResources, max_workers: int = 1) -> int:
    """
    Leave 2 cores free for the OS/other processes, then split what's left
    evenly across *max_workers* concurrent samples; cap at 24 per worker
    (diminishing returns for BWA/HaplotypeCaller beyond this in practice).

    *max_workers* matters: each ProcessPoolExecutor worker runs its own BWA/
    GATK subprocesses with this same thread count, so sizing for a single
    worker's exclusive use of the machine (the old behavior) meant
    `--max-workers 4` requested 4x the machine's actual core count. Found
    via audit.
    """
    per_worker = (resources.cpu_count - 2) // max(1, max_workers)
    return max(1, min(per_worker, 24))


def recommend_java_opts(resources: SystemResources, max_workers: int = 1) -> str:
    """
    Use 60% of *available* (not total) RAM, split evenly across
    *max_workers* concurrent JVMs so they don't collectively request more
    heap than the machine actually has (same oversubscription bug as
    `recommend_threads` — found via audit). Floors at 4g, caps at 80g per
    worker.
    """
    if resources.available_ram_gb <= 0:
        return "-Xmx8g"  # unknown RAM (non-Linux) — conservative fallback
    gb = max(4, min(int(resources.available_ram_gb * 0.6 / max(1, max_workers)), 80))
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

    Two naming conventions are accepted, per sample::

        <sample_id>_1.fastq.gz   / <sample_id>_2.fastq.gz
        <sample_id>_R1.fastq.gz  / <sample_id>_R2.fastq.gz

    Any base name matching ``*_1.fastq.gz`` or ``*_R1.fastq.gz`` is considered
    a sample; its ``_2``/``_R2`` partner is expected to exist (validated later
    per-step, see ``resolve_fastq_pair``).
    """
    pattern = re.compile(r"^(.+)_R?1\.fastq\.gz$")
    samples: set[str] = set()

    for f in input_dir.glob("*.fastq.gz"):
        m = pattern.match(f.name)
        if m:
            samples.add(m.group(1))

    if not samples:
        raise FileNotFoundError(
            f"No paired FASTQ files (matching *_1.fastq.gz or *_R1.fastq.gz) "
            f"found in {input_dir}"
        )

    return sorted(samples)


def resolve_fastq_pair(input_dir: Path, sample: str) -> tuple[Path, Path]:
    """
    Return the (R1, R2) FASTQ paths for *sample*, trying ``_1``/``_2`` first
    and falling back to ``_R1``/``_R2`` — the two conventions accepted by
    ``detect_samples``.
    """
    for suffix1, suffix2 in (("_1.fastq.gz", "_2.fastq.gz"), ("_R1.fastq.gz", "_R2.fastq.gz")):
        r1 = input_dir / f"{sample}{suffix1}"
        r2 = input_dir / f"{sample}{suffix2}"
        if r1.exists() and r2.exists():
            return r1, r2

    raise FileNotFoundError(
        f"No paired FASTQ files for sample '{sample}' in {input_dir} "
        f"(tried {sample}_1.fastq.gz/{sample}_2.fastq.gz and "
        f"{sample}_R1.fastq.gz/{sample}_R2.fastq.gz)"
    )


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
    ``<checkpoint_dir>/<sample>.<step>.<variant>.done`` exists, where
    *variant* encodes genome_build/mode/joint_genotyping.

    All three are part of the filename (not just an internal detail)
    because they change what a step's output actually *is*:
    - genome_build: switching `--genome-build` on an already-used `--output`
      directory would otherwise see the old build's checkpoints as "already
      done" and skip straight to calling variants against the new build's
      reference using a BAM still aligned to the old one.
    - joint_genotyping: HaplotypeCaller writes a different output file
      (`.g.vcf.gz` vs `.vcf`) depending on this flag; toggling it on an
      existing directory would otherwise see "haplotype" as already done
      and never produce the GVCF the cohort phase needs.
    - mode: germline and somatic filtering both write to the same
      `_PASS.vcf` filename; toggling `--mode` on an existing directory
      would otherwise leave annotation checkpointed against a PASS VCF that
      was since silently overwritten by the other mode's calls.
    Found via audit — genome_build was fixed first; mode/joint_genotyping
    had the identical gap.
    """

    def __init__(
        self,
        checkpoint_dir: Path,
        genome_build: str = "hg38",
        mode: str = "germline",
        joint_genotyping: bool = False,
    ) -> None:
        self._dir = checkpoint_dir
        self._dir.mkdir(parents=True, exist_ok=True)
        jg = "jg" if joint_genotyping else "std"
        self._variant = f"{genome_build}.{mode}.{jg}"

    def _path(self, sample: str, step: str) -> Path:
        return self._dir / f"{sample}.{step}.{self._variant}.done"

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
