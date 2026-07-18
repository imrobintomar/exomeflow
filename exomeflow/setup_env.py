"""
ExomeFlow - Automated environment setup for the ExomeFlow pipeline.

Run via:
    exomeflow setup

What happens:
  1. Auto-detects bundled GATK and ANNOVAR 
  2. Installs missing system tools (bwa, samtools, fastp, perl) 
  3. Asks the user for hg38 reference paths - or downloads them
  4. Asks the user for ANNOVAR humandb path - or downloads databases
  5. Saves everything to ~/.exomeflow/config.json for zero-arg `exomeflow run`
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import signal
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from exomeflow.config import ANNOVAR_BUILDVER
from exomeflow.utils import _parse_version, _version_ok

console = Console()
logger = logging.getLogger("exomeflow")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

GATK_VERSION = "4.6.2.0"
GATK_URL = (
    f"https://github.com/broadinstitute/gatk/releases/download/"
    f"{GATK_VERSION}/gatk-{GATK_VERSION}.zip"
)
GATK_CACHE_DIR = Path.home() / ".exomeflow" / "gatk"

# Broad Institute GATK resource bundles - publicly accessible (no auth needed)
# NOTE: hg19/GRCh37 filenames below follow the Broad legacy b37 bundle layout
# as published at implementation time; verify against the live bucket listing
# (`gsutil ls gs://gcp-public-data--broad-references/hg19/v0/`) before relying
# on them for a production GRCh37 run - flagged as an open item in the V2 plan.
_GCS_BASE_BY_BUILD: dict[str, str] = {
    "hg38":   "gs://gcp-public-data--broad-references/hg38/v0",
    "GRCh37": "gs://gcp-public-data--broad-references/hg19/v0",
}
_HTTPS_BASE_BY_BUILD: dict[str, str] = {
    "hg38":   "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0",
    "GRCh37": "https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0",
}

REFERENCE_FILES_BY_BUILD: dict[str, list[tuple[str, str, int]]] = {
    "hg38": [
        # (filename, description, approx_size_MB)
        # Core FASTA + pre-built indexes (skip `bwa index` - already in bucket)
        ("Homo_sapiens_assembly38.fasta",                       "hg38 reference genome FASTA",           3_100),
        ("Homo_sapiens_assembly38.fasta.fai",                   "hg38 FASTA index (.fai)",                   1),
        ("Homo_sapiens_assembly38.dict",                        "hg38 sequence dictionary",                  1),
        ("Homo_sapiens_assembly38.fasta.bwt",                   "BWA index (.bwt)",                        800),
        ("Homo_sapiens_assembly38.fasta.ann",                   "BWA index (.ann)",                          1),
        ("Homo_sapiens_assembly38.fasta.amb",                   "BWA index (.amb)",                          1),
        ("Homo_sapiens_assembly38.fasta.pac",                   "BWA index (.pac)",                        200),
        ("Homo_sapiens_assembly38.fasta.sa",                    "BWA index (.sa)",                         400),
        # dbSNP (correct filename in Broad bundle is dbsnp138, not dbsnp_146)
        ("Homo_sapiens_assembly38.dbsnp138.vcf.gz",             "dbSNP 138 VCF (bgzipped)",              9_500),
        ("Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi",         "dbSNP 138 tabix index",                     3),
        # BQSR known sites
        ("Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",    "Mills + 1000G gold standard indels",      130),
        ("Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi","Mills indels tabix index",                  1),
        ("Homo_sapiens_assembly38.known_indels.vcf.gz",         "Known indels VCF",                         80),
        ("Homo_sapiens_assembly38.known_indels.vcf.gz.tbi",     "Known indels tabix index",                  1),
    ],
    "GRCh37": [
        ("Homo_sapiens_assembly19.fasta",                       "GRCh37/b37 reference genome FASTA",     3_000),
        ("Homo_sapiens_assembly19.fasta.fai",                   "GRCh37 FASTA index (.fai)",                 1),
        ("Homo_sapiens_assembly19.dict",                        "GRCh37 sequence dictionary",                1),
        ("Homo_sapiens_assembly19.fasta.bwt",                   "BWA index (.bwt)",                        800),
        ("Homo_sapiens_assembly19.fasta.ann",                   "BWA index (.ann)",                          1),
        ("Homo_sapiens_assembly19.fasta.amb",                   "BWA index (.amb)",                          1),
        ("Homo_sapiens_assembly19.fasta.pac",                   "BWA index (.pac)",                        200),
        ("Homo_sapiens_assembly19.fasta.sa",                    "BWA index (.sa)",                         400),
        ("dbsnp_138.b37.vcf.gz",                                "dbSNP 138 VCF (bgzipped)",              9_000),
        ("dbsnp_138.b37.vcf.gz.tbi",                            "dbSNP 138 tabix index",                     3),
        ("Mills_and_1000G_gold_standard.indels.b37.vcf.gz",     "Mills + 1000G gold standard indels",      130),
        ("Mills_and_1000G_gold_standard.indels.b37.vcf.gz.tbi", "Mills indels tabix index",                  1),
        # Verified against the live bucket: no bgzipped "1000G_phase1.indels"
        # file exists for b37 - this is the actual known-indels file Broad
        # publishes (plain VCF + .idx, not bgzip+tbi like the other builds).
        ("Homo_sapiens_assembly19.known_indels.vcf",            "Known indels VCF",                         80),
        ("Homo_sapiens_assembly19.known_indels.vcf.idx",        "Known indels index",                        1),
    ],
}

# Mutect2's --germline-resource and --panel-of-normals - both are GATK's own
# public best-practices resources (no registration required, unlike ANNOVAR/
# OMIM). GRCh37's are dramatically larger than hg38's (a WGS-scale panel/AF
# resource vs. hg38's exome-oriented ones) - sizes verified against the live
# bucket's Content-Length at implementation time, not estimated.
_SOMATIC_BASE_BY_BUILD: dict[str, str] = {
    "hg38":   "https://storage.googleapis.com/gatk-best-practices/somatic-hg38",
    "GRCh37": "https://storage.googleapis.com/gatk-best-practices/somatic-b37",
}
_SOMATIC_GCS_BASE_BY_BUILD: dict[str, str] = {
    "hg38":   "gs://gatk-best-practices/somatic-hg38",
    "GRCh37": "gs://gatk-best-practices/somatic-b37",
}

# (filename, index_filename, description, approx_size_MB)
SOMATIC_RESOURCES_BY_BUILD: dict[str, dict[str, tuple[str, str, str, int]]] = {
    "hg38": {
        "germline_resource": (
            "af-only-gnomad.hg38.vcf.gz", "af-only-gnomad.hg38.vcf.gz.tbi",
            "gnomAD AF-only VCF (Mutect2 --germline-resource)", 3_200,
        ),
        "panel_of_normals": (
            "1000g_pon.hg38.vcf.gz", "1000g_pon.hg38.vcf.gz.tbi",
            "1000 Genomes Panel of Normals (Mutect2 --panel-of-normals)", 17,
        ),
    },
    "GRCh37": {
        "germline_resource": (
            "af-only-gnomad.raw.sites.vcf", "af-only-gnomad.raw.sites.vcf.idx",
            "gnomAD AF-only VCF (Mutect2 --germline-resource)", 14_000,
        ),
        "panel_of_normals": (
            "Mutect2-WGS-panel-b37.vcf", "Mutect2-WGS-panel-b37.vcf.idx",
            "WGS Panel of Normals (Mutect2 --panel-of-normals)", 730,
        ),
    },
}

# Alternate filenames users may already have on disk (checked for both builds)
_REF_ALTERNATES: dict[str, list[str]] = {
    "Homo_sapiens_assembly38.fasta": [
        "hg38.fa", "hg38.fasta", "GRCh38.fa", "GRCh38.fasta", "hg38.p14.fa",
    ],
    "Homo_sapiens_assembly38.dbsnp138.vcf.gz": [
        "dbsnp.vcf.gz", "dbsnp138.hg38.vcf.gz", "dbsnp_138.hg38.vcf.gz",
        "dbsnp_146.hg38.vcf.gz",   # older name some users may have
    ],
    "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz": [
        "mills.vcf.gz", "mills_indels.vcf.gz",
    ],
    "Homo_sapiens_assembly38.known_indels.vcf.gz": [
        "known_indels.vcf.gz",
    ],
    "Homo_sapiens_assembly19.fasta": [
        "hg19.fa", "hg19.fasta", "GRCh37.fa", "GRCh37.fasta", "human_g1k_v37.fasta",
    ],
    "dbsnp_138.b37.vcf.gz": [
        "dbsnp.vcf.gz", "dbsnp_138.hg19.vcf.gz",
    ],
    "Mills_and_1000G_gold_standard.indels.b37.vcf.gz": [
        "mills.vcf.gz", "mills_indels.vcf.gz",
    ],
    "Homo_sapiens_assembly19.known_indels.vcf": [
        "known_indels.vcf", "known_indels.vcf.gz",
    ],
}

# Verified against ANNOVAR's own `-webfrom annovar` mirror listing
# (hg38_avdblist.txt / hg19_avdblist.txt) at implementation time - every
# entry here is confirmed actually downloadable through
# `_download_annovar_db`'s mechanism, and matches config.py's
# ANNOVAR_PROTOCOLS_BY_BUILD exactly (db names must stay in sync).
# COSMIC is deliberately excluded: it requires separate registered download
# and was never obtainable through this mirror.
ANNOVAR_DATABASES_BY_BUILD: dict[str, list[tuple[str, str, int]]] = {
    "hg38": [
        # (db_name, description, approx_size_MB)
        ("refGene",            "Gene annotation (RefSeq)",               80),
        ("avsnp150",           "dbSNP 150 with allele frequencies",  10_000),
        ("clinvar_20240611",   "ClinVar clinical significance",          120),
        ("gnomad41_exome",     "gnomAD v4.1 exome allele frequencies", 8_000),
        ("gnomad41_genome",    "gnomAD v4.1 genome allele frequencies",30_000),
        ("dbnsfp47a",          "dbNSFP 4.7a functional predictions",  50_000),
        ("exac03",             "ExAC 0.3 allele frequencies",          1_200),
    ],
    "GRCh37": [
        # gnomAD v4.1 was never released for hg19/GRCh37 - v2.1.1 is the
        # newest gnomAD build available there. Everything else matches hg38.
        ("refGene",            "Gene annotation (RefSeq)",               80),
        ("avsnp150",           "dbSNP 150 with allele frequencies",  10_000),
        ("clinvar_20240611",   "ClinVar clinical significance",          120),
        ("gnomad211_exome",    "gnomAD v2.1.1 exome allele frequencies", 500),
        ("gnomad211_genome",   "gnomAD v2.1.1 genome allele frequencies",5_500),
        ("dbnsfp47a",          "dbNSFP 4.7a functional predictions",  50_000),
        ("exac03",             "ExAC 0.3 allele frequencies",          1_200),
    ],
}

_SYSTEM_TOOLS: list[tuple[str, str]] = [
    # (binary_name, conda_spec)
    ("bwa",      "bioconda::bwa"),
    ("samtools", "bioconda::samtools"),
    ("fastp",    "bioconda::fastp"),
    ("perl",     "conda-forge::perl"),
]

# micromamba - a single self-contained binary, no installer/shell-init
# required - bootstrapped automatically when neither conda nor mamba is on
# PATH, so bwa/samtools/fastp/perl can still be auto-installed without
# requiring the user to go install Miniconda themselves first.
MICROMAMBA_DIR = Path.home() / ".exomeflow" / "conda"
MICROMAMBA_BIN = MICROMAMBA_DIR / "bin" / "micromamba"
MICROMAMBA_ENV = MICROMAMBA_DIR / "env"
MICROMAMBA_URL = "https://micro.mamba.pm/api/micromamba/linux-64/latest"

CONFIG_PATH = Path.home() / ".exomeflow" / "config.json"

# Canonical (reference, dbsnp, mills, known_indels) filenames per build -
# must match the first/dbSNP/Mills/known-indels entries in REFERENCE_FILES_BY_BUILD.
_CANONICAL_REF_NAMES: dict[str, dict[str, str]] = {
    "hg38": {
        "reference": "Homo_sapiens_assembly38.fasta",
        "dbsnp": "Homo_sapiens_assembly38.dbsnp138.vcf.gz",
        "mills": "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        "known_indels": "Homo_sapiens_assembly38.known_indels.vcf.gz",
    },
    "GRCh37": {
        "reference": "Homo_sapiens_assembly19.fasta",
        "dbsnp": "dbsnp_138.b37.vcf.gz",
        "mills": "Mills_and_1000G_gold_standard.indels.b37.vcf.gz",
        "known_indels": "Homo_sapiens_assembly19.known_indels.vcf",
    },
}

# ---------------------------------------------------------------------------
# Config file
# ---------------------------------------------------------------------------

def load_config() -> dict:
    """Load saved ExomeFlow config. Returns empty dict if not found."""
    if CONFIG_PATH.exists():
        try:
            return json.loads(CONFIG_PATH.read_text())
        except Exception:
            return {}
    return {}


def save_config(data: dict) -> None:
    """Merge data into the saved ExomeFlow config."""
    existing = load_config()
    existing.update({k: str(v) for k, v in data.items() if v is not None})
    CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    CONFIG_PATH.write_text(json.dumps(existing, indent=2))


# ---------------------------------------------------------------------------
# Auto-detection
# ---------------------------------------------------------------------------

def detect_gatk_bin() -> Path | None:
    """
    Return path to the 'gatk' executable.
    Checks the saved config first (so a previously-resolved path survives
    running the CLI from a different working directory), then ExomeFlow
    source root, PATH, and common locations.
    """
    saved = load_config().get("gatk_bin")
    if saved and Path(saved).is_file():
        return Path(saved)

    source_root = Path(__file__).parent.parent
    for name in [f"gatk-{GATK_VERSION}", "gatk"]:
        p = source_root / name / "gatk"
        if p.is_file():
            return p

    found = shutil.which("gatk")
    if found:
        return Path(found)

    for candidate in [
        Path.home() / "gatk" / "gatk",
        Path.home() / f"gatk-{GATK_VERSION}" / "gatk",
        Path("/opt/gatk/gatk"),
        Path(f"/opt/gatk-{GATK_VERSION}/gatk"),
        GATK_CACHE_DIR / f"gatk-{GATK_VERSION}" / "gatk",
    ]:
        if candidate.is_file():
            return candidate

    return None


def _step_gatk_download() -> Path | None:
    """
    Auto-download and extract GATK into ~/.exomeflow/gatk/ if not found
    anywhere. Needs Java on PATH to actually run afterward - GATK's zip
    bundles jars + a wrapper script, not a JVM, which pip cannot provision.
    """
    logger.info("Downloading GATK %s (~600 MB) ...", GATK_VERSION)
    GATK_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    zip_path = GATK_CACHE_DIR / f"gatk-{GATK_VERSION}.zip"

    if not _download_file(GATK_URL, zip_path):
        logger.warning("GATK download failed - install manually: %s", GATK_URL)
        return None

    try:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(GATK_CACHE_DIR)
    except Exception as exc:
        logger.warning("GATK archive extraction failed: %s", exc)
        return None
    finally:
        zip_path.unlink(missing_ok=True)

    gatk_bin = GATK_CACHE_DIR / f"gatk-{GATK_VERSION}" / "gatk"
    if not gatk_bin.is_file():
        logger.warning("GATK extracted but 'gatk' script not found at %s", gatk_bin)
        return None
    gatk_bin.chmod(gatk_bin.stat().st_mode | 0o111)

    if not shutil.which("java"):
        logger.warning(
            "GATK downloaded, but no 'java' on PATH - GATK needs Java 17+ to "
            "run. Install a JDK (e.g. `conda install -c conda-forge openjdk=17`)."
        )

    logger.log(25, "GATK installed: %s", gatk_bin)
    return gatk_bin


def detect_annovar_bin() -> Path | None:
    """
    Return path to the ANNOVAR directory containing table_annovar.pl.
    Checks the saved config first (so a previously-resolved path survives
    running the CLI from a different working directory), then ExomeFlow
    source root, then common locations.
    """
    saved = load_config().get("annovar_bin")
    if saved and (Path(saved) / "table_annovar.pl").exists():
        return Path(saved)

    candidates = [
        Path(__file__).parent.parent / "annovar",
        Path.cwd() / "annovar",
        Path.home() / "annovar",
        Path("/opt/annovar"),
        Path("/usr/local/annovar"),
    ]
    for p in candidates:
        if (p / "table_annovar.pl").exists():
            return p
    return None


def detect_annovar_humandb(buildver: str = "hg38", timeout_s: int = 30) -> Path | None:
    """
    Look for an existing ANNOVAR humandb directory anywhere on the system -
    databases are 10s-100s of GB, so reusing one already downloaded (e.g. by
    another tool, or a prior ExomeFlow install) beats re-downloading blindly.

    Fast path: a short list of common locations. Fallback: a bounded `find`
    across mounted filesystems (maxdepth-limited, so it stays fast even on
    multi-TB drives) - skipped entirely if `find` isn't available.
    """
    fast_candidates = [
        Path.home() / ".exomeflow" / "annovar" / "humandb",
        Path(__file__).parent.parent / "annovar" / "humandb",
        Path.cwd() / "annovar" / "humandb",
        Path("/opt/annovar/humandb"),
        Path("/usr/local/annovar/humandb"),
    ]
    for p in fast_candidates:
        if (p / f"{buildver}_refGene.txt").exists():
            return p

    if not shutil.which("find"):
        return None
    roots = [str(p) for p in (Path("/media"), Path("/mnt"), Path("/data"), Path.home()) if p.is_dir()]
    if not roots:
        return None
    stdout = ""
    try:
        result = subprocess.run(
            ["find", *roots, "-maxdepth", "6", "-iname", f"{buildver}_refGene.txt"],
            capture_output=True, text=True, timeout=timeout_s,
        )
        stdout = result.stdout
    except subprocess.TimeoutExpired as exc:
        # Use whatever find had already printed before the timeout instead
        # of discarding a real hit just because the rest of the scan (e.g.
        # one more, slower mount) didn't finish in time.
        stdout = exc.stdout or ""

    hits = [Path(line).parent for line in stdout.splitlines() if line]
    if not hits:
        return None
    if len(hits) == 1:
        return hits[0]

    # More than one refGene.txt on the system is possible - e.g. InterVar
    # bundles its own small humandb subset that also has a refGene file -
    # and `find`'s output order isn't meaningful. Prefer whichever hit has
    # the most {buildver}_*.txt files alongside it, as a proxy for "the
    # real, full humandb" over a small bundled subset.
    def _file_count(p: Path) -> int:
        try:
            return sum(1 for _ in p.glob(f"{buildver}_*.txt"))
        except OSError:
            return 0

    return max(hits, key=_file_count)


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def _run_bounded(cmd: list[str], timeout_s: int, *, visible: bool) -> bool:
    """
    Run *cmd* with a hard wall-clock timeout, killing the whole process
    group (not just the direct child) if it's exceeded.

    Found via audit: a plain subprocess.run(timeout=...) only kills the
    direct child - several commands here (perl's annotate_variation.pl,
    git) can themselves spawn a grandchild download process that survives
    the parent's death and keeps running/consuming bandwidth forever. The
    same class of bug was already found and fixed for InterVar
    (acmg_classification.py); this generalizes that fix to every other
    subprocess call in this module that was still unguarded.
    """
    kwargs = {} if visible else {"stdout": subprocess.DEVNULL, "stderr": subprocess.DEVNULL}
    proc = subprocess.Popen(cmd, text=True, start_new_session=True, **kwargs)
    try:
        proc.wait(timeout=timeout_s)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()
        logger.warning("Command timed out after %ds and was killed: %s", timeout_s, " ".join(cmd))
        return False
    return proc.returncode == 0


def _run_silent(cmd: list[str], timeout_s: int = 900) -> bool:
    return _run_bounded(cmd, timeout_s, visible=False)


def _run_visible(cmd: list[str], timeout_s: int = 21600) -> bool:
    return _run_bounded(cmd, timeout_s, visible=True)


def _ask(prompt: str, default_yes: bool = False, assume_yes: bool = False) -> bool:
    if assume_yes or not sys.stdin.isatty():
        # --yes, or no interactive terminal to prompt on at all (background/CI
        # runs) - fall back to the default instead of hanging on input().
        return True if assume_yes else default_yes
    suffix = "[Y/n]" if default_yes else "[y/N]"
    try:
        ans = input(f"  {prompt} {suffix}: ").strip().lower()
        if not ans:
            return default_yes
        return ans in ("y", "yes")
    except (EOFError, KeyboardInterrupt):
        return default_yes


def _ask_path(prompt: str) -> Path | None:
    try:
        raw = input(f"  {prompt}: ").strip()
        if raw:
            p = Path(raw).expanduser()
            return p if p.exists() else None
        return None
    except (EOFError, KeyboardInterrupt):
        return None


def _download_file(url: str, dest: Path) -> bool:
    """
    Download with wget (resume-capable) or curl. These are legitimately
    multi-hour operations for a ~3GB reference FASTA, so no fixed outer
    timeout is right here - instead, use each tool's own stall-detection
    (abort if the *transfer rate*, not the whole operation, stalls), plus a
    generous 6-hour outer safety net via _run_visible's default.
    """
    for tool, cmd in [
        ("wget",  ["wget", "-c", "-q", "--show-progress", "--timeout=60", "--tries=3",
                    url, "-O", str(dest)]),
        ("curl",  ["curl", "-L", "--progress-bar", "-C", "-",
                    "--speed-time", "60", "--speed-limit", "1024",
                    url, "-o", str(dest)]),
    ]:
        if shutil.which(tool):
            return _run_visible(cmd)
    console.print("  [red]✘[/red]  Neither wget nor curl found.")
    return False


def _gsutil_cp(gcs_path: str, dest: Path) -> bool:
    return _run_visible(["gsutil", "-m", "cp", gcs_path, str(dest)])


# ---------------------------------------------------------------------------
# Step 1 — Bundled tools (GATK + ANNOVAR)
# ---------------------------------------------------------------------------

def _step_bundled_tools(assume_yes: bool = False) -> tuple[Path | None, Path | None]:
    """
    Detect bundled GATK and ANNOVAR. Add GATK to PATH.
    Returns (gatk_path, annovar_bin).
    """
    console.print(Panel("[bold]Step 1 - Bundled Tools[/bold]", style="blue"))

    # GATK - auto-download into ~/.exomeflow/gatk/ if not found anywhere
    gatk = detect_gatk_bin() or _step_gatk_download()
    if gatk:
        try:
            gatk.chmod(gatk.stat().st_mode | 0o111)
        except Exception:
            pass
        gatk_dir = str(gatk.parent)
        if gatk_dir not in os.environ.get("PATH", ""):
            os.environ["PATH"] = gatk_dir + os.pathsep + os.environ.get("PATH", "")
        console.print(f"  [green]✔[/green]  GATK {GATK_VERSION}: [cyan]{gatk}[/cyan]")
    else:
        console.print(
            f"  [red]✘[/red]  GATK auto-download failed.\n"
            f"  Place gatk-{GATK_VERSION}/ inside the ExomeFlow folder or download from:\n"
            f"  [cyan]{GATK_URL}[/cyan]"
        )

    # ANNOVAR - cannot be auto-downloaded: it requires free personal
    # registration before Annovar's own site will hand out a download link
    # (their license, not ours). Once you have it, it can live anywhere.
    annovar = detect_annovar_bin()
    if annovar:
        console.print(f"  [green]✔[/green]  ANNOVAR: [cyan]{annovar}[/cyan]")
    else:
        console.print(
            "  [red]✘[/red]  ANNOVAR not found automatically.\n"
            "  Register (free) and download at: "
            "[cyan]https://annovar.openbioinformatics.org[/cyan]\n"
            "  Extract the tar.gz anywhere on disk - it doesn't need to be inside\n"
            "  any particular folder."
        )
        if not assume_yes and _ask("Already have ANNOVAR downloaded and extracted somewhere?"):
            user_path = _ask_path("Enter the path to the directory containing table_annovar.pl")
            if user_path and (user_path / "table_annovar.pl").exists():
                annovar = user_path
                console.print(f"  [green]✔[/green]  ANNOVAR: [cyan]{annovar}[/cyan]")
            else:
                console.print("  [red]✘[/red]  table_annovar.pl not found at that path.")

    return gatk, annovar


# ---------------------------------------------------------------------------
# Step 2 — System tools
# ---------------------------------------------------------------------------

def _step_bootstrap_micromamba() -> Path | None:
    """
    Download a self-contained micromamba binary into ~/.exomeflow/conda/ 
    used only when neither conda nor mamba is already on PATH, so
    bwa/samtools/fastp/perl can still be auto-installed without requiring
    the user to install Miniconda themselves first. No installer, no shell
    init: it's a single ~7 MB static binary.
    """
    if MICROMAMBA_BIN.is_file():
        return MICROMAMBA_BIN

    console.print("  [yellow]→[/yellow]  conda/mamba not found - bootstrapping micromamba (~7 MB) ...")
    MICROMAMBA_BIN.parent.mkdir(parents=True, exist_ok=True)
    archive = MICROMAMBA_DIR / "micromamba.tar.bz2"

    if not _download_file(MICROMAMBA_URL, archive):
        console.print("  [red]✘[/red]  micromamba download failed.")
        return None

    try:
        with tarfile.open(archive) as tf:
            member = tf.getmember("bin/micromamba")
            member.name = MICROMAMBA_BIN.name
            tf.extract(member, MICROMAMBA_BIN.parent)
    except Exception as exc:
        console.print(f"  [red]✘[/red]  micromamba archive extraction failed: {exc}")
        return None
    finally:
        archive.unlink(missing_ok=True)

    if not MICROMAMBA_BIN.is_file():
        console.print(f"  [red]✘[/red]  micromamba extracted but binary not found at {MICROMAMBA_BIN}")
        return None

    MICROMAMBA_BIN.chmod(MICROMAMBA_BIN.stat().st_mode | 0o111)
    console.print(f"  [green]✔[/green]  micromamba installed: [cyan]{MICROMAMBA_BIN}[/cyan]")
    return MICROMAMBA_BIN


def _step_system_tools(outdated: frozenset[str] = frozenset()) -> list[str]:
    """
    Install missing system tools via conda/mamba, or a self-bootstrapped
    micromamba if neither is present; force-upgrade any already on PATH but
    below the required minimum version (*outdated*). Returns list of
    failures.
    """
    console.print(Panel("[bold]Step 2 - System Tools[/bold]", style="blue"))
    failures = []

    conda = shutil.which("conda") or shutil.which("mamba")
    using_micromamba = False
    if not conda:
        bootstrapped = _step_bootstrap_micromamba()
        if not bootstrapped:
            console.print(
                "  [red]✘[/red]  Could not auto-install a package manager.\n"
                "  Install Miniconda manually: [cyan]https://docs.conda.io/en/latest/miniconda.html[/cyan]"
            )
            return ["conda not found and micromamba bootstrap failed - install tools manually"]
        conda = str(bootstrapped)
        using_micromamba = True
        MICROMAMBA_ENV.mkdir(parents=True, exist_ok=True)
        env_bin = str(MICROMAMBA_ENV / "bin")
        if env_bin not in os.environ.get("PATH", ""):
            os.environ["PATH"] = env_bin + os.pathsep + os.environ.get("PATH", "")

    for binary, spec in _SYSTEM_TOOLS:
        if shutil.which(binary) and binary not in outdated:
            console.print(f"  [green]✔[/green]  {binary} already installed")
            continue
        verb = "Upgrading" if binary in outdated else "Installing"
        console.print(f"  [yellow]→[/yellow]  {verb} {binary} ...")
        channel, pkg = spec.split("::")
        if using_micromamba:
            # "create" (not "install") - micromamba's "install" refuses to
            # target a prefix that isn't already a created environment, but
            # "create" against an existing prefix just adds the package
            # alongside what's already there (verified: safe to call
            # repeatedly against the same fixed MICROMAMBA_ENV prefix).
            cmd = [conda, "create", "-y", "-p", str(MICROMAMBA_ENV),
                   "-c", channel, "-c", "conda-forge", pkg]
        else:
            cmd = [conda, "install", "-y", "-c", channel, pkg]
        ok = _run_silent(cmd, timeout_s=1800)
        if ok:
            console.print(f"  [green]✔[/green]  {binary} installed")
        else:
            console.print(f"  [red]✘[/red]  {binary} install failed - try: conda install -c {channel} {pkg}")
            failures.append(f"{binary} install failed")

    return failures


# ---------------------------------------------------------------------------
# Step 3 - Reference files
# ---------------------------------------------------------------------------

def _scan_refs(directory: Path, genome_build: str = "hg38") -> dict[str, Path]:
    """Return {canonical_name: found_path} for all reference files found in directory."""
    found: dict[str, Path] = {}
    for filename, _, _ in REFERENCE_FILES_BY_BUILD[genome_build]:
        p = directory / filename
        if p.exists():
            found[filename] = p
            continue
        for alt in _REF_ALTERNATES.get(filename, []):
            p = directory / alt
            if p.exists():
                found[filename] = p
                break
    return found


def _step_reference_files(
    refs_dir: Path, existing_refs_dir: Path | None, genome_build: str = "hg38",
    assume_yes: bool = False,
) -> tuple[Path | None, list[str]]:
    """
    Locate or download reference files for *genome_build*.
    Returns (resolved_refs_dir, failures).
    """
    reference_files = REFERENCE_FILES_BY_BUILD[genome_build]
    required_names = set(_CANONICAL_REF_NAMES[genome_build].values())
    console.print(Panel(f"[bold]Step 3 - Reference Genome ({genome_build})[/bold]", style="blue"))
    failures: list[str] = []

    # If user passed --existing-refs, check there first. Requires all 4
    # actually-needed files (fasta + dbSNP/Mills/known-indels VCFs) - not
    # just *any* single match among the 14 possible files/indexes
    # _scan_refs looks for. Found via audit: a directory holding only a
    # FASTA (no VCFs at all - a very plausible real scenario for a user who
    # already has a reference genome but not GATK's known-sites bundle) used
    # to be accepted as complete, so the required VCFs were never fetched.
    if existing_refs_dir:
        found = _scan_refs(existing_refs_dir, genome_build)
        if required_names.issubset(found):
            console.print(f"  [green]✔[/green]  Found {len(found)}/{len(reference_files)} reference files in [cyan]{existing_refs_dir}[/cyan]")
            return existing_refs_dir, []
        missing = required_names - set(found)
        console.print(
            f"  [yellow]⚠[/yellow]  {existing_refs_dir} is missing required file(s): "
            f"{', '.join(sorted(missing))} — will resolve the rest separately."
        )

    # Check default location - same all-4-required check as above.
    found = _scan_refs(refs_dir, genome_build)
    if required_names.issubset(found):
        console.print(f"  [green]✔[/green]  Reference files already present in [cyan]{refs_dir}[/cyan]")
        return refs_dir, []

    # Ask user
    console.print(
        f"  {genome_build} reference files are required (~13 GB total).\n"
        "  These include: reference FASTA, dbSNP, Mills indels, known indels."
    )

    # --yes never answers this one: it needs a manually-entered path, so with
    # no path to give it, "no, go download instead" is the only sane default.
    if not assume_yes and _ask(f"Do you already have {genome_build} reference files on this machine?"):
        user_path = _ask_path("Enter the path to the directory containing your reference files")
        if user_path:
            found = _scan_refs(user_path, genome_build)
            if required_names.issubset(found):
                console.print(f"  [green]✔[/green]  Found {len(found)} reference file(s) in [cyan]{user_path}[/cyan]")
                return user_path, []
            elif found:
                # Partial match: same all-4-required rule as the other two
                # checks above (found via audit). Rather than discarding a
                # path the user explicitly typed in, fill in the missing
                # piece(s) at that same location instead of falling back to
                # the unrelated default refs_dir.
                missing = required_names - set(found)
                console.print(
                    f"  [yellow]⚠[/yellow]  Found {len(found)} reference file(s) in [cyan]{user_path}[/cyan], "
                    f"but missing: {', '.join(sorted(missing))} — will download the rest there."
                )
                refs_dir = user_path
            else:
                console.print(f"  [red]✘[/red]  No recognised reference files found in {user_path}")
                failures.append(f"No reference files found in {user_path}")
                return None, failures
        else:
            console.print("  [red]✘[/red]  Path not found or not entered.")
            failures.append("Reference path not provided")
            return None, failures

    # Offer download
    total_mb = sum(s for _, _, s in reference_files)
    if not _ask(
        f"Download reference files to {refs_dir}? (~{total_mb // 1024} GB, may take hours)",
        default_yes=True, assume_yes=assume_yes,
    ):
        failures.append("Reference files not downloaded - required for pipeline")
        return None, failures

    refs_dir.mkdir(parents=True, exist_ok=True)
    use_gsutil = bool(shutil.which("gsutil"))
    if not use_gsutil:
        console.print(
            "  [dim]gsutil not found - using wget. "
            "For faster downloads: conda install -c conda-forge google-cloud-sdk[/dim]"
        )

    gcs_base = _GCS_BASE_BY_BUILD[genome_build]
    https_base = _HTTPS_BASE_BY_BUILD[genome_build]

    for filename, description, size_mb in reference_files:
        dest = refs_dir / filename
        if dest.exists() and dest.stat().st_size > 0:
            console.print(f"  [green]✔[/green]  {filename} already present")
            continue
        console.print(f"  [cyan]→[/cyan]  Downloading {filename} ({description}, ~{size_mb:,} MB) ...")
        # Not trusting the tool's own return value alone: gsutil's sliced/
        # resumable-download bookkeeping can report a spurious non-zero
        # exit even when the file actually transferred completely and
        # correctly (found live, fixed for somatic resources in 2.2.7 —
        # applying the same fix here, since this download loop had the
        # identical pattern).
        if use_gsutil:
            _gsutil_cp(f"{gcs_base}/{filename}", dest)
        else:
            _download_file(f"{https_base}/{filename}", dest)
        if dest.exists() and dest.stat().st_size > 0:
            console.print(f"  [green]✔[/green]  {filename}")
        else:
            console.print(f"  [red]✘[/red]  {filename} download failed")
            failures.append(f"Download failed: {filename}")

    return refs_dir if not failures else None, failures


def _step_somatic_resources(
    refs_dir: Path, genome_build: str = "hg38", assume_yes: bool = False,
) -> dict[str, Path]:
    """
    Auto-resolve/download Mutect2's --germline-resource and
    --panel-of-normals for --mode somatic. Both are GATK's own public
    best-practices resources — no registration required, unlike ANNOVAR/OMIM.

    Best-effort and never raises: Mutect2 still runs without either (with a
    higher false-positive rate, per run_mutect2()'s own warning) — these are
    accuracy accelerants, not hard requirements, so a failed/declined
    download here just means the run proceeds without that resource rather
    than blocking.
    """
    resources = SOMATIC_RESOURCES_BY_BUILD[genome_build]
    https_base = _SOMATIC_BASE_BY_BUILD[genome_build]
    gcs_base = _SOMATIC_GCS_BASE_BY_BUILD[genome_build]
    use_gsutil = bool(shutil.which("gsutil"))

    resolved: dict[str, Path] = {}
    for key, (filename, index_filename, description, size_mb) in resources.items():
        dest = refs_dir / filename
        dest_idx = refs_dir / index_filename
        if dest.exists() and dest.stat().st_size > 0 and dest_idx.exists() and dest_idx.stat().st_size > 0:
            resolved[key] = dest
            continue

        # The germline-resource AF file is multi-GB (14 GB for GRCh37) — ask
        # before pulling that much data. The PoN (17-730 MB) just downloads,
        # matching the HPO mapping's small-file auto-download pattern.
        if size_mb > 1_000 and not _ask(
            f"Download {description} (~{size_mb // 1024 or 1} GB)?",
            default_yes=True, assume_yes=assume_yes,
        ):
            logger.info("Skipping %s - Mutect2 will run without it.", description)
            continue

        refs_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Downloading %s (~%s MB) ...", description, f"{size_mb:,}")
        # Deliberately not `and`-chained: gsutil's sliced/resumable-download
        # bookkeeping can report a spurious non-zero exit (e.g. failing to
        # rename a temp component file) even when the actual transfer
        # completed correctly - found live, where the main VCF downloaded
        # byte-perfect but a false failure signal on it meant the .tbi
        # index was never even attempted. Both are always tried
        # independently; success is judged by the files actually existing
        # and being non-empty afterward, not by trusting the tool's exit
        # code for either.
        if use_gsutil:
            _gsutil_cp(f"{gcs_base}/{filename}", dest)
            _gsutil_cp(f"{gcs_base}/{index_filename}", dest_idx)
        else:
            _download_file(f"{https_base}/{filename}", dest)
            _download_file(f"{https_base}/{index_filename}", dest_idx)
        ok = (
            dest.exists() and dest.stat().st_size > 0
            and dest_idx.exists() and dest_idx.stat().st_size > 0
        )
        if ok:
            logger.log(25, "%s ready.", description)
            resolved[key] = dest
        else:
            logger.warning("%s download failed - Mutect2 will run without it.", description)

    return resolved


# ---------------------------------------------------------------------------
# Step 4 - ANNOVAR databases
# ---------------------------------------------------------------------------

def annovar_databases_complete(annovar_db: Path, genome_build: str = "hg38") -> tuple[bool, list[str]]:
    """
    Check whether every database `cfg.annovar_protocols` actually needs is
    present in *annovar_db* for *genome_build* - not just that the directory
    exists. Returns (complete, missing_db_names).

    Restored via audit: the steady-state pre-flight check had narrowed to a
    bare directory-exists check, so a database deleted after initial setup
    (or a build mismatch) went undetected until table_annovar.pl failed
    hours into a run instead of being caught in seconds up front.

    Filter-type ("f" operation) databases also need their paired .idx file 
    ANNOVAR's `-downdb -webfrom annovar` fetches both files together, but a
    connection drop partway through the paired download can land the .txt
    without its .idx. refGene is gene-based ("g" operation) and doesn't use
    this index format, so it's exempt.
    """
    buildver = ANNOVAR_BUILDVER[genome_build]
    required = ANNOVAR_DATABASES_BY_BUILD[genome_build]
    missing = []
    for d, _, _ in required:
        txt = annovar_db / f"{buildver}_{d}.txt"
        if not txt.exists():
            missing.append(d)
        elif d != "refGene" and not Path(str(txt) + ".idx").exists():
            missing.append(d)
    return not missing, missing


def _step_annovar_databases(
    annovar_bin: Path, default_db: Path, genome_build: str = "hg38",
    assume_yes: bool = False,
) -> tuple[Path | None, list[str]]:
    """
    Locate or download ANNOVAR databases for *genome_build* ("hg38" or "GRCh37").
    Returns (humandb_path, failures).
    """
    console.print(Panel("[bold]Step 4 - ANNOVAR Annotation Databases[/bold]", style="blue"))
    failures: list[str] = []
    buildver = ANNOVAR_BUILDVER[genome_build]
    required = ANNOVAR_DATABASES_BY_BUILD[genome_build]

    if not annovar_bin:
        console.print("  [yellow]⚠[/yellow]  ANNOVAR not detected  skipping database step.")
        failures.append("ANNOVAR not found  databases not set up")
        return None, failures

    annotate_pl = annovar_bin / "annotate_variation.pl"
    total_gb = sum(s for _, _, s in required) // 1024

    def _report(label: str, path: Path, missing: list[str], auto: bool = False) -> None:
        found_n = len(required) - len(missing)
        suffix = " (auto-detected)" if auto else ""
        if not missing:
            console.print(
                f"  [green]✔[/green]  Found {found_n}/{len(required)} "
                f"databases in [cyan]{path}[/cyan]{suffix}"
            )
        else:
            console.print(
                f"  [yellow]⚠[/yellow]  Found {found_n}/{len(required)} "
                f"databases in [cyan]{path}[/cyan]{suffix} — missing: {', '.join(missing)}"
            )

    # Check existing db in the annovar_bin/humandb location. A location that
    # has *some* but not all required databases isn't returned as done — it
    # falls through to the download loop below so only the missing ones get
    # fetched. (Found via a live run: a partial match here used to be
    # reported as complete and the pipeline proceeded to annotate without
    # ever downloading the database it had just told the user was missing.)
    db_dir: Path | None = None
    # dict.fromkeys dedupes while preserving order — default_db and
    # annovar_bin/humandb are frequently the same path (whenever no
    # explicit --annovar-db was given), which used to print the identical
    # "Found N/7 ... missing: ..." line twice in a row.
    for candidate_db in dict.fromkeys([default_db, annovar_bin / "humandb"]):
        if candidate_db.exists():
            complete, missing = annovar_databases_complete(candidate_db, genome_build)
            _report("default", candidate_db, missing)
            if complete:
                return candidate_db, []
            if len(missing) < len(required):
                db_dir = candidate_db
                break

    if db_dir is None:
        # System-wide lookup these are 10s-100s of GB; check for an existing
        # humandb elsewhere on disk before asking the user or downloading blind.
        console.print("  [dim]Searching the system for an existing ANNOVAR humandb ...[/dim]")
        found_db = detect_annovar_humandb(buildver)
        if found_db:
            complete, missing = annovar_databases_complete(found_db, genome_build)
            _report("system", found_db, missing, auto=True)
            if complete:
                return found_db, []
            db_dir = found_db

    if db_dir is None:
        # Ask user
        console.print(
            f"  ANNOVAR databases are required for variant annotation (~{total_gb} GB total).\n"
            "  These include: refGene, ClinVar, gnomAD, dbNSFP, avSNP."
        )

        if not assume_yes and _ask("Do you already have an ANNOVAR humandb directory?"):
            user_path = _ask_path("Enter the path to your humandb directory")
            if user_path and user_path.is_dir():
                complete, missing = annovar_databases_complete(user_path, genome_build)
                _report("user-supplied", user_path, missing)
                if complete:
                    return user_path, []
                db_dir = user_path
            else:
                console.print("  [red]✘[/red]  Path not found or not entered.")
                failures.append("ANNOVAR humandb path not valid")
                return None, failures

    # Offer download (either a fresh directory, or filling in the gaps in a
    # partially-complete one found above)
    if db_dir is None:
        db_dir = default_db
    if not _ask(
        f"Download the missing ANNOVAR database(s) to {db_dir}? (~{total_gb} GB total set)",
        default_yes=True, assume_yes=assume_yes,
    ):
        failures.append("ANNOVAR databases not downloaded — required for annotation")
        return None, failures

    db_dir.mkdir(parents=True, exist_ok=True)

    for db_name, description, size_mb in required:
        db_file = db_dir / f"{buildver}_{db_name}.txt"
        # refGene is gene-based ("g" operation) and has no .idx pair; every
        # other database here is filter-type ("f") and needs both files 
        # re-download (ANNOVAR fetches the pair together) if either is
        # missing, matching annovar_databases_complete()'s own check.
        idx_ok = db_name == "refGene" or Path(str(db_file) + ".idx").exists()
        if db_file.exists() and idx_ok:
            console.print(f"  [green]✔[/green]  {db_name} already present")
            continue
        console.print(f"  [cyan]→[/cyan]  Downloading {db_name} ({description}, ~{size_mb:,} MB) ...")
        # 3h timeout: dbnsfp47a alone is ~50GB. Guards the same class of
        # orphaned-grandchild-process bug already fixed for InterVar, since
        # annotate_variation.pl -downdb can itself spawn a sub-downloader.
        _run_visible(
            ["perl", str(annotate_pl), "-buildver", buildver,
             "-downdb", "-webfrom", "annovar", db_name, str(db_dir)],
            timeout_s=10800,
        )
        # Not trusting the exit code alone (same class of bug already fixed
        # for gsutil elsewhere in this module)  verify the actual files
        # landed and are non-empty.
        idx_file = Path(str(db_file) + ".idx")
        idx_ok = db_name == "refGene" or (idx_file.exists() and idx_file.stat().st_size > 0)
        if db_file.exists() and db_file.stat().st_size > 0 and idx_ok:
            console.print(f"  [green]✔[/green]  {db_name}")
        else:
            console.print(f"  [red]✘[/red]  {db_name} download failed")
            failures.append(f"ANNOVAR db download failed: {db_name}")

    return db_dir if not failures else None, failures


# ---------------------------------------------------------------------------
# Step 5  InterVar (ACMG classification) + HPO gene-to-phenotype mapping
# ---------------------------------------------------------------------------

INTERVAR_URL = "https://github.com/WGLab/InterVar/archive/refs/heads/master.zip"
INTERVAR_DIR = Path.home() / ".exomeflow" / "intervar"


def detect_intervar_bin() -> Path | None:
    """Return the directory containing Intervar.py  bundled folder first, then cache."""
    candidates = [
        Path(__file__).parent.parent / "intervar",
        Path(__file__).parent.parent / "InterVar",
        Path.cwd() / "intervar",
        INTERVAR_DIR,
        INTERVAR_DIR / "InterVar-master",
    ]
    for p in candidates:
        if (p / "Intervar.py").exists():
            return p
    return None


MIM2GENE_URL = "https://omim.org/static/omim/data/mim2gene.txt"


def _ensure_mim2gene(intervar_dir: Path) -> None:
    """
    InterVar's ACMG classification step hard-requires intervardb/mim2gene.txt
    to run at all  without it, InterVar prints "Error: can't read the OMIM
    file ... Please download it from http://www.omim.org/downloads" and
    silently fails to produce any classification output (the pipeline then
    only sees a generic "expected output not found" with no indication why).

    Unlike ANNOVAR, this specific file needs no registration  OMIM
    publishes it as a plain public download  but InterVar's own GitHub repo
    doesn't ship it (likely the same redistribution concern that stops
    ANNOVAR bundling ClinVar), so a git clone alone never provisions it.
    Found via a live run: ACMG classification silently failed on every
    sample despite InterVar itself being correctly installed.
    """
    dest = intervar_dir / "intervardb" / "mim2gene.txt"
    if dest.exists():
        return
    logger.info("Downloading OMIM mim2gene.txt (required by InterVar's ACMG classification) ...")
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not _download_file(MIM2GENE_URL, dest):
        logger.warning(
            "mim2gene.txt download failed  ACMG classification will fail until "
            "it's placed at %s (freely downloadable, no registration required).",
            dest,
        )


def _step_intervar(annovar_bin: Path | None) -> Path | None:
    """Auto-clone InterVar and download its database if not already present."""
    found = detect_intervar_bin()
    if found:
        _ensure_mim2gene(found)
        return found

    if not shutil.which("git"):
        logger.warning("git not found  cannot auto-install InterVar for ACMG classification.")
        return None

    logger.info("Installing InterVar (ACMG classification) into %s ...", INTERVAR_DIR)
    INTERVAR_DIR.parent.mkdir(parents=True, exist_ok=True)
    ok = _run_silent(
        ["git", "clone", "--depth", "1",
         "https://github.com/WGLab/InterVar.git", str(INTERVAR_DIR)]
    )
    if not ok or not (INTERVAR_DIR / "Intervar.py").exists():
        logger.warning("InterVar auto-install failed  ACMG classification will be skipped.")
        return None

    _ensure_mim2gene(INTERVAR_DIR)
    logger.log(25, "InterVar installed: %s", INTERVAR_DIR)
    return INTERVAR_DIR


def _step_hpo_mapping() -> bool:
    """Auto-download the HPO gene-to-phenotype mapping used for HPO enrichment."""
    from exomeflow.hpo_annotation import HPO_CACHE_DIR, HPO_DOWNLOAD_URL, HPO_MAPPING_FILE

    if HPO_MAPPING_FILE.exists():
        return True
    HPO_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading HPO gene-to-phenotype mapping ...")
    ok = _download_file(HPO_DOWNLOAD_URL, HPO_MAPPING_FILE)
    if not ok:
        logger.warning("HPO mapping download failed  HPO enrichment will be skipped.")
    return ok


def _step_multiqc() -> bool:
    """Auto-pip-install multiqc for the cohort rollup report."""
    if shutil.which("multiqc"):
        return True
    logger.info("Installing multiqc via pip ...")
    ok = _run_silent([sys.executable, "-m", "pip", "install", "multiqc"])
    if not ok:
        logger.warning("multiqc auto-install failed  MultiQC rollup will be skipped.")
        return False
    return True


def _step_matplotlib() -> bool:
    """Auto-pip-install matplotlib, only needed when --cnv is requested."""
    try:
        import matplotlib  # noqa: F401
        return True
    except ImportError:
        pass
    logger.info("Installing matplotlib via pip (required for --cnv plots) ...")
    ok = _run_silent([sys.executable, "-m", "pip", "install", "matplotlib>=3.7.0"])
    if not ok:
        logger.warning("matplotlib auto-install failed  CNV plotting may fail.")
        return False
    return True


# ---------------------------------------------------------------------------
# Pre-flight dependency check (called on every `exomeflow run`)
# ---------------------------------------------------------------------------

_TOOL_CHECKS: list[tuple[str, list[str], str, str]] = [
    # (name, version_cmd, label, min_version) — min_version restored via
    # audit: this table used to only carry a description, and the check
    # loop below had narrowed to presence-only, silently dropping the
    # minimum-version enforcement the old (now-deleted) requirements
    # checker used to do.
    ("bwa",      ["bwa"],                    "BWA aligner",      "0.7.17"),
    ("samtools", ["samtools", "--version"],   "SAMtools",         "1.13"),
    ("fastp",    ["fastp", "--version"],      "fastp QC",         "0.20.1"),
    ("perl",     ["perl", "--version"],       "Perl interpreter", "5.26"),
    ("gatk",     ["gatk", "--version"],       "GATK 4",           "4.6.0"),
]


def _tool_version_ok(version_cmd: list[str], min_version: str) -> tuple[bool, str]:
    """Run *version_cmd* and check its output against *min_version*. Returns (ok, found_version_str)."""
    try:
        result = subprocess.run(
            version_cmd, capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=10,
        )
    except (subprocess.TimeoutExpired, OSError):
        return False, "?"
    found = _parse_version(result.stdout + result.stderr)
    return _version_ok(found, min_version), ".".join(str(v) for v in found)

_REF_KEYS: list[tuple[str, str]] = [
    # (config_key, label)
    ("reference",    "Reference FASTA"),
    ("dbsnp",        "dbSNP VCF"),
    ("mills",        "Mills indels VCF"),
    ("known_indels", "Known indels VCF"),
]

_ANNOVAR_KEYS: list[tuple[str, str]] = [
    ("annovar_bin", "ANNOVAR scripts directory"),
    ("annovar_db",  "ANNOVAR humandb directory"),
]


def _ask_genome_build(saved_build: str | None, assume_yes: bool) -> str:
    """Resolve which genome build to use: saved config > explicit prompt > hg38 default."""
    if saved_build in ("hg38", "GRCh37"):
        return saved_build
    if assume_yes or not sys.stdin.isatty():
        return "hg38"
    try:
        ans = input("  Which genome build? [hg38/GRCh37] (default hg38): ").strip()
    except (EOFError, KeyboardInterrupt):
        return "hg38"
    return "GRCh37" if ans.lower() in ("grch37", "hg19", "37") else "hg38"


def check_and_fix_dependencies(
    genome_build: str | None = None, call_cnv: bool = False, assume_yes: bool = False,
    mode: str = "germline",
) -> dict:
    """
    Run a pre-flight check before every pipeline run.

    Prints a status table showing every dependency.
    If anything is missing: runs the relevant fix step automatically.
    Returns the final config dict (with all paths resolved) on success.
    Raises SystemExit(1) if any required *hard* dependency cannot be fixed.
    MultiQC / HPO / InterVar are best-effort and never raise — see Phase 3/7
    of the V2 plan: their absence degrades output, it doesn't block the run.

    *genome_build*: None means "not explicitly set by the user on the CLI" —
    resolved from saved config if present, otherwise prompted for (or
    defaulted to hg38 under --yes/non-interactive).
    """
    from rich.panel import Panel

    # Python version floor — dropped silently when the old requirements
    # checker was deleted; restored here. Bug found via audit.
    if sys.version_info[:2] < (3, 9):
        py_ver = ".".join(str(v) for v in sys.version_info[:3])
        console.print(
            f"  [red]✘[/red]  Python {py_ver} < 3.9 — ExomeFlow requires Python "
            f"3.9 or newer. Upgrade: conda install python>=3.9"
        )
        raise SystemExit(1)

    cfg = load_config()
    saved_build = cfg.get("genome_build")
    genome_build = _ask_genome_build(genome_build or saved_build, assume_yes)
    cfg["genome_build"] = genome_build
    # Switching build on an already-configured install: the saved reference/
    # ANNOVAR paths point at the *old* build's files, which still exist on
    # disk (a plain .exists() check can't tell), so a naive check would
    # silently keep running the old build's refs against the new build's
    # buildver — force re-resolution instead. Bug found via audit.
    build_switched = bool(saved_build) and saved_build != genome_build

    # ── Build status for all checks ─────────────────────────────────────────
    issues: list[str] = []   # categories that need fixing

    table = Table(show_header=True, header_style="bold cyan", box=None, padding=(0, 2))
    table.add_column("Category",    style="bold white",  min_width=22)
    table.add_column("Dependency",  style="white",       min_width=30)
    table.add_column("Status",      justify="center",    min_width=12)

    # Tools
    # GATK resolved once up front (bundled-folder detection counts as
    # "present" even when not literally on PATH yet) so the printed table
    # and the issues list agree — previously the table always showed GATK
    # as "missing" whenever it wasn't on PATH, then a later block silently
    # patched `issues` without the already-printed table ever reflecting
    # it, contradicting what the code had actually decided. Found via audit.
    gatk_path = detect_gatk_bin()
    tool_missing: list[str] = []
    tool_outdated: list[str] = []
    for name, version_cmd, label, min_version in _TOOL_CHECKS:
        present = bool(shutil.which(name)) or (name == "gatk" and gatk_path is not None)
        version_ok, found_ver = (True, "")
        if present and shutil.which(name):
            # Only version-check a tool actually resolved via PATH — a
            # bundled-but-not-yet-PATH'd GATK is checked once it's added to
            # PATH later in this same function.
            version_ok, found_ver = _tool_version_ok(version_cmd, min_version)
        ok = present and version_ok
        if present and shutil.which(name) and not version_ok:
            status = f"[yellow]⚠ {found_ver} < {min_version}[/yellow]"
            tool_outdated.append(name)
        else:
            status = "[green]✔ found[/green]" if ok else "[red]✘ missing[/red]"
        table.add_row("Tools" if name == "bwa" else "", f"{label} ({name})", status)
        if not ok:
            tool_missing.append(name)
    if tool_missing:
        issues.append("tools")

    # Reference files
    ref_missing: list[str] = []
    for key, label in _REF_KEYS:
        val = cfg.get(key)
        ok = bool(val and Path(val).exists()) and not build_switched
        status = "[green]✔ found[/green]" if ok else "[red]✘ missing[/red]"
        table.add_row("References" if key == "reference" else "", label, status)
        if not ok:
            ref_missing.append(key)
    if ref_missing:
        issues.append("refs")

    # ANNOVAR
    annovar_missing: list[str] = []
    for key, label in _ANNOVAR_KEYS:
        val = cfg.get(key)
        # annovar_bin itself isn't build-specific (only annovar_db's contents
        # are), so build_switched only forces re-resolution of annovar_db.
        ok = bool(val and Path(val).exists()) and not (build_switched and key == "annovar_db")
        # Restored via audit: a directory existing doesn't mean every database
        # this build's protocol list needs is actually inside it — check the
        # specific files too, so a database deleted after initial setup (or
        # simply never fetched for this build) is caught here in seconds
        # instead of failing table_annovar.pl hours into a real run.
        if ok and key == "annovar_db":
            complete, missing_dbs = annovar_databases_complete(Path(val), genome_build)
            if not complete:
                ok = False
                console.print(
                    f"  [yellow]⚠[/yellow]  ANNOVAR humandb missing database(s): "
                    f"{', '.join(missing_dbs)}"
                )
        status = "[green]✔ found[/green]" if ok else "[red]✘ missing[/red]"
        table.add_row("ANNOVAR" if key == "annovar_bin" else "", label, status)
        if not ok:
            annovar_missing.append(key)
    if annovar_missing:
        issues.append("annovar")

    if build_switched:
        console.print(
            f"  [yellow]⚠  Genome build changed ({saved_build} → {genome_build}) — "
            "re-resolving reference/ANNOVAR paths for the new build.[/yellow]"
        )

    # ── Print the table ──────────────────────────────────────────────────────
    console.print()
    console.print(Panel(table, title="[bold]ExomeFlow — Dependency Check[/bold]",
                        border_style="blue", expand=False))

    if not issues:
        console.print("  [green]✔  All dependencies satisfied.[/green] Starting pipeline...\n")
        if saved_build != genome_build:
            save_config({"genome_build": genome_build})  # persist even on the fast path
        _step_multiqc()
        _step_hpo_mapping()
        _step_intervar(detect_annovar_bin())
        if call_cnv:
            _step_matplotlib()
        if mode == "somatic":
            refs_dir = Path(cfg.get("refs_dir", Path.home() / ".exomeflow" / "refs"))
            somatic_paths = _step_somatic_resources(refs_dir, genome_build, assume_yes)
            if somatic_paths:
                save_config({k: str(v) for k, v in somatic_paths.items()})
                cfg.update({k: str(v) for k, v in somatic_paths.items()})
        return cfg

    # ── Fix what's missing ──────────────────────────────────────────────────
    console.print(f"\n  [yellow]⚠  Fixing {len(issues)} missing requirement(s)...[/yellow]\n")

    new_cfg: dict = {}

    # Fix tools
    if "tools" in issues:
        console.print(Panel("[bold]Fixing: System Tools[/bold]", style="yellow"))
        # GATK — add bundled path to PATH, or auto-download if not found anywhere
        if "gatk" in tool_missing:
            resolved_gatk = gatk_path or _step_gatk_download()
            if resolved_gatk:
                _add_to_path(resolved_gatk)
                console.print(f"  [green]✔[/green]  GATK resolved: [cyan]{resolved_gatk}[/cyan]")
                tool_missing.remove("gatk")
                new_cfg["gatk_bin"] = resolved_gatk
            else:
                console.print(
                    f"  [red]✘[/red]  GATK auto-download failed. Place "
                    f"[cyan]gatk-{GATK_VERSION}/[/cyan] inside the ExomeFlow folder, "
                    f"or download manually: {GATK_URL}"
                )
        # Other tools via conda (installs anything missing, upgrades anything outdated)
        remaining_tools = [t for t in tool_missing if t != "gatk"]
        if remaining_tools:
            _step_system_tools(outdated=frozenset(tool_outdated))

        # A tool that's still missing, or still below its minimum version
        # after an upgrade attempt, is a hard failure — previously this fell
        # through silently and the pipeline died later on a raw "command not
        # found" (or an obscure tool-internal error for outdated versions)
        # instead of a clear upfront error.
        still_missing = [t for t in tool_missing if not shutil.which(t)]
        for name, version_cmd, _label, min_version in _TOOL_CHECKS:
            if name in tool_outdated and shutil.which(name):
                ok, found_ver = _tool_version_ok(version_cmd, min_version)
                if not ok:
                    still_missing.append(f"{name} {found_ver} < {min_version} (upgrade failed)")
        if still_missing:
            console.print(f"  [red]✘[/red]  Still missing: {', '.join(still_missing)}")
            raise SystemExit(1)

    # Fix reference files
    if "refs" in issues:
        console.print(Panel("[bold]Fixing: Reference Files[/bold]", style="yellow"))
        refs_dir = Path(cfg.get("refs_dir", Path.home() / ".exomeflow" / "refs"))
        resolved_refs, ref_failures = _step_reference_files(
            refs_dir, existing_refs_dir=None, genome_build=genome_build,
            assume_yes=assume_yes,
        )
        if resolved_refs is None:
            console.print("  [red]✘[/red]  Reference files are required. Cannot continue.")
            raise SystemExit(1)
        names = _CANONICAL_REF_NAMES[genome_build]
        new_cfg["refs_dir"]     = resolved_refs
        new_cfg["reference"]    = _find_ref(resolved_refs, names["reference"])
        new_cfg["dbsnp"]        = _find_ref(resolved_refs, names["dbsnp"])
        new_cfg["mills"]        = _find_ref(resolved_refs, names["mills"])
        new_cfg["known_indels"] = _find_ref(resolved_refs, names["known_indels"])

    # Fix ANNOVAR
    if "annovar" in issues:
        console.print(Panel("[bold]Fixing: ANNOVAR[/bold]", style="yellow"))
        annovar_bin = detect_annovar_bin()
        if annovar_bin is None:
            console.print(
                "  [red]✘[/red]  ANNOVAR not found automatically.\n"
                "  Register (free) and download at: "
                "[cyan]https://annovar.openbioinformatics.org[/cyan]\n"
                "  Extract the tar.gz anywhere on disk — it doesn't need to be inside\n"
                "  any particular folder."
            )
            if not assume_yes and _ask("Already have ANNOVAR downloaded and extracted somewhere?"):
                user_path = _ask_path("Enter the path to the directory containing table_annovar.pl")
                if user_path and (user_path / "table_annovar.pl").exists():
                    annovar_bin = user_path
                    console.print(f"  [green]✔[/green]  ANNOVAR: [cyan]{annovar_bin}[/cyan]")
        if annovar_bin is None:
            console.print("  [red]✘[/red]  ANNOVAR is required. Cannot continue.")
            raise SystemExit(1)
        default_db = Path(cfg.get("annovar_db", annovar_bin / "humandb"))
        resolved_db, db_failures = _step_annovar_databases(
            annovar_bin, default_db, genome_build, assume_yes=assume_yes
        )
        if resolved_db is None:
            console.print("  [red]✘[/red]  ANNOVAR databases are required. Cannot continue.")
            raise SystemExit(1)
        new_cfg["annovar_bin"] = annovar_bin
        new_cfg["annovar_db"]  = resolved_db

    # Save any newly resolved paths
    new_cfg["genome_build"] = genome_build
    if new_cfg:
        save_config(new_cfg)
        cfg = load_config()
        console.print(f"\n  [dim]Configuration updated → {CONFIG_PATH}[/dim]")

    # Best-effort extras — never block the run (Phase 3/7 of the V2 plan)
    _step_multiqc()
    _step_hpo_mapping()
    _step_intervar(detect_annovar_bin())
    if call_cnv:
        _step_matplotlib()
    if mode == "somatic":
        refs_dir = Path(cfg.get("refs_dir", Path.home() / ".exomeflow" / "refs"))
        somatic_paths = _step_somatic_resources(refs_dir, genome_build, assume_yes)
        if somatic_paths:
            save_config({k: str(v) for k, v in somatic_paths.items()})
            cfg.update({k: str(v) for k, v in somatic_paths.items()})

    console.print(Panel("[bold green]✔  All dependencies satisfied. Starting pipeline...[/bold green]",
                        style="green"))
    console.print()
    return cfg


def _add_to_path(executable: Path) -> None:
    """Make executable's parent directory the first entry on PATH."""
    try:
        executable.chmod(executable.stat().st_mode | 0o111)
    except Exception:
        pass
    d = str(executable.parent)
    if d not in os.environ.get("PATH", ""):
        os.environ["PATH"] = d + os.pathsep + os.environ.get("PATH", "")


def _find_ref(directory: Path, canonical: str) -> Path:
    """Return the first existing alternate for a canonical reference filename."""
    for name in [canonical] + _REF_ALTERNATES.get(canonical, []):
        p = directory / name
        if p.exists():
            return p
    return directory / canonical


# ---------------------------------------------------------------------------
# Read-only pre-flight report (`exomeflow doctor`) — never downloads,
# installs, or writes config; safe to run any time.
# ---------------------------------------------------------------------------

def run_doctor(genome_build: str | None = None) -> None:
    """
    Print one consolidated report of every dependency ExomeFlow needs:
    found/missing, and — critically — whether a missing one will resolve
    itself automatically on the next `exomeflow run`/`exomeflow setup`, or
    needs the user to go do something first (this only ever applies to
    ANNOVAR, which requires free registration before it can be downloaded —
    every other dependency here auto-resolves).

    Purely diagnostic: makes no network requests, installs nothing, and
    never writes ~/.exomeflow/config.json.
    """
    from rich.panel import Panel

    cfg = load_config()
    genome_build = genome_build or cfg.get("genome_build") or "hg38"

    table = Table(show_header=True, header_style="bold cyan", box=None, padding=(0, 2))
    table.add_column("Category",      style="bold white", min_width=22)
    table.add_column("Dependency",    style="white",       min_width=30)
    table.add_column("Status",        justify="center",    min_width=14)
    table.add_column("If missing",    style="dim",          min_width=28)

    py_ver = ".".join(str(v) for v in sys.version_info[:3])
    py_ok = sys.version_info[:2] >= (3, 9)
    table.add_row(
        "Python", f"Interpreter ({py_ver})",
        "[green]✔ found[/green]" if py_ok else "[red]✘ too old[/red]",
        "" if py_ok else "Upgrade to Python >= 3.9",
    )

    conda = shutil.which("conda") or shutil.which("mamba")
    conda_status = "[green]✔ found[/green]" if conda else (
        "[yellow]⚠ bootstrapped[/yellow]" if MICROMAMBA_BIN.is_file() else "[yellow]⚠ missing[/yellow]"
    )
    table.add_row(
        "Package manager", "conda / mamba (for bwa/samtools/fastp/perl)", conda_status,
        "" if conda or MICROMAMBA_BIN.is_file()
        else "Auto-bootstraps a local micromamba — no action needed",
    )

    gatk_path = detect_gatk_bin()
    for name, version_cmd, label, min_version in _TOOL_CHECKS:
        if name == "gatk":
            present, found_ver = gatk_path is not None, ""
        else:
            present = bool(shutil.which(name))
            found_ver = ""
        version_ok = True
        if present and (shutil.which(name) or name == "gatk"):
            check_cmd = version_cmd if shutil.which(name) else None
            if check_cmd:
                version_ok, found_ver = _tool_version_ok(check_cmd, min_version)
        ok = present and version_ok
        if present and found_ver and not version_ok:
            status = f"[yellow]⚠ {found_ver} < {min_version}[/yellow]"
        else:
            status = "[green]✔ found[/green]" if ok else "[red]✘ missing[/red]"
        table.add_row(
            "Tools" if name == "bwa" else "", f"{label} ({name})", status,
            "" if ok else "Auto-installed on next setup/run",
        )

    for key, label in _REF_KEYS:
        val = cfg.get(key)
        ok = bool(val and Path(val).exists())
        table.add_row(
            "References" if key == "reference" else "", label,
            "[green]✔ found[/green]" if ok else "[red]✘ missing[/red]",
            "" if ok else "Auto-downloaded on next setup/run",
        )

    annovar_bin = detect_annovar_bin()
    table.add_row(
        "ANNOVAR", "Scripts directory (annovar_bin)",
        "[green]✔ found[/green]" if annovar_bin else "[red]✘ missing[/red]",
        "" if annovar_bin else "Register + download yourself — see below",
    )
    # detect_annovar_humandb() expects ANNOVAR's own buildver naming
    # ("hg38"/"hg19"), not the genome_build the rest of the CLI uses
    # ("hg38"/"GRCh37") — every other call site converts via
    # ANNOVAR_BUILDVER first. Found via audit: passing the raw genome_build
    # through here meant a GRCh37 user with no annovar_db saved yet always
    # searched for "GRCh37_refGene.txt" (which never exists) instead of
    # "hg19_refGene.txt", so `exomeflow doctor` reported ANNOVAR databases
    # as missing even when a fully-populated humandb was on disk.
    annovar_db = Path(cfg["annovar_db"]) if cfg.get("annovar_db") else (
        detect_annovar_humandb(ANNOVAR_BUILDVER[genome_build]) if annovar_bin else None
    )
    db_ok = False
    if annovar_db and annovar_db.exists():
        db_ok, _ = annovar_databases_complete(annovar_db, genome_build)
    table.add_row(
        "", "Annotation databases (annovar_db)",
        "[green]✔ found[/green]" if db_ok else "[red]✘ missing[/red]",
        "" if db_ok else ("Auto-downloaded once ANNOVAR itself is set up" if annovar_bin
                           else "Needs ANNOVAR itself first"),
    )

    intervar = detect_intervar_bin()
    table.add_row(
        "Optional extras", "InterVar (ACMG classification)",
        "[green]✔ found[/green]" if intervar else "[yellow]⚠ missing[/yellow]",
        "" if intervar else "Auto-installed; degrades gracefully if it can't be",
    )
    from exomeflow.hpo_annotation import HPO_MAPPING_FILE
    hpo_ok = HPO_MAPPING_FILE.exists()
    table.add_row(
        "", "HPO gene-phenotype mapping",
        "[green]✔ found[/green]" if hpo_ok else "[yellow]⚠ missing[/yellow]",
        "" if hpo_ok else "Auto-downloaded; degrades gracefully if it can't be",
    )
    multiqc_ok = bool(shutil.which("multiqc"))
    table.add_row(
        "", "MultiQC (cohort QC rollup)",
        "[green]✔ found[/green]" if multiqc_ok else "[yellow]⚠ missing[/yellow]",
        "" if multiqc_ok else "Auto-installed via pip; degrades gracefully if it can't be",
    )

    console.print()
    console.print(Panel(table, title="[bold]ExomeFlow — Doctor[/bold]",
                        border_style="blue", expand=False))

    if not annovar_bin:
        console.print(
            "\n  [yellow]⚠  ANNOVAR is the one dependency that can't be auto-installed[/yellow] "
            "— it requires free\n"
            "  personal registration before its own site will hand out a download link "
            "(that's\n"
            "  ANNOVAR's license, not a limitation of ExomeFlow). Everything else in this "
            "report\n"
            "  resolves itself automatically the next time you run "
            "[cyan]exomeflow setup[/cyan] or\n"
            "  [cyan]exomeflow run[/cyan].\n\n"
            "  1. Register + download: [cyan]https://annovar.openbioinformatics.org[/cyan]\n"
            "  2. Extract the tar.gz anywhere on disk\n"
            "  3. Run [cyan]exomeflow setup[/cyan] — it will ask for the path "
            "interactively\n"
        )
    else:
        console.print(
            "\n  Run [cyan]exomeflow setup[/cyan] to resolve anything still missing above, "
            "or just\n  [cyan]exomeflow run[/cyan] — first-run setup happens automatically.\n"
        )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_setup(
    refs_dir: Path,
    annovar_bin: Path | None = None,
    annovar_db: Path | None = None,
    existing_refs_dir: Path | None = None,
    genome_build: str | None = None,
    assume_yes: bool = False,
) -> int:
    """
    Run the full setup sequence.
    Returns number of failures (0 = all good).
    """
    console.print(Panel(
        "[bold white]ExomeFlow — Setup[/bold white]\n"
        "[dim]Installs tools, configures references and annotation databases.[/dim]",
        style="bold blue",
        expand=False,
    ))

    all_failures: list[str] = []
    config: dict = {}
    saved_cfg = load_config()
    genome_build = _ask_genome_build(genome_build or saved_cfg.get("genome_build"), assume_yes)

    # ── Step 1: Bundled GATK + ANNOVAR ──────────────────────────────────────
    gatk_path, resolved_annovar_bin = _step_bundled_tools(assume_yes=assume_yes)
    if gatk_path:
        config["gatk_bin"] = gatk_path
    effective_annovar_bin = annovar_bin or resolved_annovar_bin

    # ── Step 2: System tools ────────────────────────────────────────────────
    all_failures += _step_system_tools()

    # ── Step 3: Reference files ─────────────────────────────────────────────
    resolved_refs, ref_failures = _step_reference_files(
        refs_dir, existing_refs_dir, genome_build=genome_build, assume_yes=assume_yes
    )
    all_failures += ref_failures
    if resolved_refs:
        names = _CANONICAL_REF_NAMES[genome_build]
        config["refs_dir"]     = resolved_refs
        config["genome_build"] = genome_build
        config["reference"]    = _find_ref(resolved_refs, names["reference"])
        config["dbsnp"]        = _find_ref(resolved_refs, names["dbsnp"])
        config["mills"]        = _find_ref(resolved_refs, names["mills"])
        config["known_indels"] = _find_ref(resolved_refs, names["known_indels"])

    # ── Step 4: ANNOVAR databases ────────────────────────────────────────────
    if effective_annovar_bin:
        # A previously-saved annovar_db (e.g. answered interactively in an
        # earlier `exomeflow setup`/`exomeflow run`) is checked before
        # falling back to <annovar_bin>/humandb — otherwise every re-run of
        # `exomeflow setup` "forgot" it and re-triggered the whole
        # find-or-ask flow from scratch. Not trusted blindly: still run
        # through _step_annovar_databases()'s own existence/completeness
        # check like any other candidate.
        saved_db = saved_cfg.get("annovar_db")
        default_db = annovar_db or (Path(saved_db) if saved_db else effective_annovar_bin / "humandb")
        resolved_db, db_failures = _step_annovar_databases(
            effective_annovar_bin, default_db, genome_build, assume_yes=assume_yes
        )
        all_failures += db_failures
        if resolved_db:
            config["annovar_bin"] = effective_annovar_bin
            config["annovar_db"]  = resolved_db

    # ── Step 5: InterVar + HPO mapping (best-effort) ─────────────────────────
    _step_intervar(effective_annovar_bin)
    _step_hpo_mapping()
    _step_multiqc()

    # ── Save config ──────────────────────────────────────────────────────────
    if config:
        save_config(config)
        console.print(f"\n  [dim]Configuration saved to {CONFIG_PATH}[/dim]")

    # ── Summary ──────────────────────────────────────────────────────────────
    console.print()
    if not all_failures:
        console.print(Panel(
            "[bold green]✔  Setup complete![/bold green]\n\n"
            "Run your pipeline:\n"
            "[cyan]exomeflow run --input-dir fastq/ --output results/[/cyan]",
            style="green",
        ))
    else:
        console.print(Panel(
            f"[bold yellow]Setup finished with {len(all_failures)} issue(s):[/bold yellow]\n\n"
            + "\n".join(f"  [red]•[/red] {f}" for f in all_failures),
            style="yellow",
        ))

    return len(all_failures)
