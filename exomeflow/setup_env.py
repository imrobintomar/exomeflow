"""
ExomeFlow — Automated environment setup.

Installs all missing tools, Python packages, ANNOVAR databases, and
human genome reference files in one command.

Run via:
    exomeflow setup --refs-dir /path/to/refs --annovar-bin /path/to/annovar --annovar-db /path/to/humandb
"""

from __future__ import annotations

import hashlib
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path

from rich.console import Console
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
    TransferSpeedColumn,
)
from rich.table import Table

logger = logging.getLogger("exomeflow")
console = Console()

# ---------------------------------------------------------------------------
# Python packages
# ---------------------------------------------------------------------------

_PYTHON_PACKAGES: list[tuple[str, str]] = [
    ("typer", "0.12.0"),
    ("rich", "13.0.0"),
    ("pandas", "2.0.0"),
    ("matplotlib", "3.7.0"),
]

# ---------------------------------------------------------------------------
# System tools
# ---------------------------------------------------------------------------

_SYSTEM_TOOLS: list[tuple[str, str, str]] = [
    # (name, conda_pkg, min_version)
    ("bwa", "bioconda::bwa", "0.7.17"),
    ("samtools", "bioconda::samtools", "1.13"),
    ("fastp", "bioconda::fastp", "0.20.1"),
    ("perl", "conda-forge::perl", "5.26"),
]

# ---------------------------------------------------------------------------
# GATK download
# ---------------------------------------------------------------------------

GATK_VERSION = "4.6.2.0"
GATK_URL = (
    f"https://github.com/broadinstitute/gatk/releases/download/"
    f"{GATK_VERSION}/gatk-{GATK_VERSION}.zip"
)

# ---------------------------------------------------------------------------
# Reference files — GATK resource bundle (GCS public bucket, no auth needed)
# ---------------------------------------------------------------------------

_GCS_BASE = "gs://genomics-public-data/references/hg38/v0"
_GATK_BASE = "gs://gatk-best-practices/somatic-hg38"

REFERENCE_FILES: list[tuple[str, str, str, int]] = [
    # (filename, gcs_path, description, approx_size_MB)
    (
        "Homo_sapiens_assembly38.fasta",
        f"{_GCS_BASE}/Homo_sapiens_assembly38.fasta",
        "hg38 reference genome FASTA",
        3_100,
    ),
    (
        "Homo_sapiens_assembly38.fasta.fai",
        f"{_GCS_BASE}/Homo_sapiens_assembly38.fasta.fai",
        "hg38 FASTA index (.fai)",
        1,
    ),
    (
        "Homo_sapiens_assembly38.dict",
        f"{_GCS_BASE}/Homo_sapiens_assembly38.dict",
        "hg38 sequence dictionary",
        1,
    ),
    (
        "dbsnp_146.hg38.vcf.gz",
        f"{_GCS_BASE}/dbsnp_146.hg38.vcf.gz",
        "dbSNP 146 VCF (bgzipped)",
        9_500,
    ),
    (
        "dbsnp_146.hg38.vcf.gz.tbi",
        f"{_GCS_BASE}/dbsnp_146.hg38.vcf.gz.tbi",
        "dbSNP tabix index",
        3,
    ),
    (
        "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        f"{_GCS_BASE}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        "Mills + 1000G gold standard indels",
        130,
    ),
    (
        "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
        f"{_GCS_BASE}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
        "Mills indels tabix index",
        1,
    ),
    (
        "Homo_sapiens_assembly38.known_indels.vcf.gz",
        f"{_GCS_BASE}/Homo_sapiens_assembly38.known_indels.vcf.gz",
        "Known indels VCF",
        80,
    ),
    (
        "Homo_sapiens_assembly38.known_indels.vcf.gz.tbi",
        f"{_GCS_BASE}/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi",
        "Known indels tabix index",
        1,
    ),
]

# ---------------------------------------------------------------------------
# ANNOVAR databases
# ---------------------------------------------------------------------------

ANNOVAR_DATABASES: list[tuple[str, str, int]] = [
    # (db_name, description, approx_size_MB)
    ("refGene", "Gene annotation (RefSeq)", 80),
    ("avsnp150", "dbSNP 150 with allele frequencies", 10_000),
    ("clinvar_20240416", "ClinVar clinical significance", 120),
    ("gnomad41_exome", "gnomAD v4.1 exome allele frequencies", 8_000),
    ("gnomad41_genome", "gnomAD v4.1 genome allele frequencies", 30_000),
    ("dbnsfp47a", "dbNSFP 4.7a functional predictions", 50_000),
    ("cosmic84_coding", "COSMIC v84 somatic mutations", 500),
    ("exac03", "ExAC 0.3 allele frequencies", 1_200),
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _ask(question: str, default_yes: bool = False) -> bool:
    """Prompt the user and return True for yes."""
    suffix = "[Y/n]" if default_yes else "[y/N]"
    try:
        answer = input(f"\n  {question} {suffix}: ").strip().lower()
        if not answer:
            return default_yes
        return answer in ("y", "yes")
    except (EOFError, KeyboardInterrupt):
        return False


def _run(cmd: list[str], desc: str = "") -> bool:
    """Run a shell command, stream output, return True on success."""
    console.print(f"  [cyan]→[/cyan] {desc or ' '.join(cmd)}")
    result = subprocess.run(cmd, text=True)
    return result.returncode == 0


def _pip_install(package: str, min_ver: str) -> bool:
    spec = f"{package}>={min_ver}"
    console.print(f"  [cyan]→[/cyan] pip install {spec}")
    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", "-q", spec],
        text=True,
    )
    return result.returncode == 0


def _conda_install(conda_pkg: str) -> bool:
    conda = shutil.which("conda") or shutil.which("mamba")
    if not conda:
        console.print("  [red]✘[/red]  conda/mamba not found.")
        return False
    channel, pkg = conda_pkg.split("::")
    result = subprocess.run(
        [conda, "install", "-y", "-c", channel, pkg],
        text=True,
    )
    return result.returncode == 0


def _gsutil_available() -> bool:
    return shutil.which("gsutil") is not None


def _download_gsutil(gcs_path: str, dest: Path) -> bool:
    """Download a single file from GCS using gsutil."""
    result = subprocess.run(
        ["gsutil", "-m", "cp", gcs_path, str(dest)],
        text=True,
    )
    return result.returncode == 0


def _download_wget(url: str, dest: Path) -> bool:
    """Download a file using wget with progress bar."""
    result = subprocess.run(
        ["wget", "-c", "-q", "--show-progress", url, "-O", str(dest)],
        text=True,
    )
    return result.returncode == 0


def _bwa_index(fasta: Path) -> bool:
    """Index the reference genome with BWA."""
    console.print("  [cyan]→[/cyan] Indexing reference with BWA (this takes ~1 hour) ...")
    result = subprocess.run(["bwa", "index", str(fasta)], text=True)
    return result.returncode == 0


def _samtools_faidx(fasta: Path) -> bool:
    """Create .fai index with samtools."""
    result = subprocess.run(["samtools", "faidx", str(fasta)], text=True)
    return result.returncode == 0


def _gatk_dict(fasta: Path) -> bool:
    """Create sequence dictionary with GATK."""
    dict_file = fasta.with_suffix(".dict")
    result = subprocess.run(
        ["gatk", "CreateSequenceDictionary", "-R", str(fasta), "-O", str(dict_file)],
        text=True,
    )
    return result.returncode == 0


# ---------------------------------------------------------------------------
# Setup steps
# ---------------------------------------------------------------------------

def setup_python_packages() -> list[str]:
    """Install / upgrade all Python dependencies. Returns list of failures."""
    console.print(Panel("[bold]Step 1 — Python Packages[/bold]", style="blue"))
    failures = []
    for pkg, min_ver in _PYTHON_PACKAGES:
        try:
            import importlib
            mod = importlib.import_module(pkg)
            ver = getattr(mod, "__version__", "?")
            console.print(f"  [green]✔[/green]  {pkg} {ver} already installed")
        except ImportError:
            console.print(f"  [yellow]⚠[/yellow]  {pkg} missing — installing ...")
            if not _pip_install(pkg, min_ver):
                failures.append(f"pip install {pkg}>={min_ver} failed")
    return failures


def setup_system_tools() -> list[str]:
    """Install missing system tools via conda. Returns list of failures."""
    console.print(Panel("[bold]Step 2 — System Tools[/bold]", style="blue"))
    failures = []

    for name, conda_pkg, min_ver in _SYSTEM_TOOLS:
        if shutil.which(name):
            console.print(f"  [green]✔[/green]  {name} already on PATH")
            continue

        console.print(f"  [yellow]⚠[/yellow]  {name} not found.")
        if _ask(f"Install {name} via conda? ({conda_pkg})"):
            if not _conda_install(conda_pkg):
                failures.append(f"conda install {conda_pkg} failed")
            else:
                console.print(
                    f"  [green]✔[/green]  {name} installed. "
                    f"[dim]Open a new terminal if it's still not on PATH.[/dim]"
                )
        else:
            failures.append(f"{name} not installed — required for pipeline")

    # GATK — cannot conda install, needs manual setup
    if not shutil.which("gatk"):
        console.print("  [yellow]⚠[/yellow]  GATK not found on PATH.")
        console.print(
            f"  [dim]Download GATK {GATK_VERSION} from:[/dim]\n"
            f"  [link]{GATK_URL}[/link]\n"
            f"  [dim]Then add to PATH:[/dim]\n"
            f"  [cyan]export PATH=/path/to/gatk-{GATK_VERSION}:$PATH[/cyan]"
        )
        failures.append("GATK must be installed manually — see instructions above")
    else:
        console.print("  [green]✔[/green]  GATK already on PATH")

    return failures


def _find_existing_refs(refs_dir: Path) -> dict[str, Path]:
    """
    Scan refs_dir for known reference filenames.
    Also accepts common alternate names (e.g. hg38.fa, hg38.fasta).
    Returns {canonical_filename: found_path}.
    """
    # Alternate names users might already have
    ALTERNATES: dict[str, list[str]] = {
        "Homo_sapiens_assembly38.fasta": [
            "hg38.fa", "hg38.fasta", "GRCh38.fa", "GRCh38.fasta",
            "hg38.p14.fa", "Homo_sapiens_assembly38.fa",
        ],
        "dbsnp_146.hg38.vcf.gz": [
            "dbsnp.vcf.gz", "dbsnp138.hg38.vcf.gz", "dbsnp_138.hg38.vcf.gz",
        ],
        "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz": [
            "mills.vcf.gz", "mills_indels.vcf.gz",
        ],
        "Homo_sapiens_assembly38.known_indels.vcf.gz": [
            "known_indels.vcf.gz", "Homo_sapiens_assembly38.known_indels.vcf.gz",
        ],
    }
    found: dict[str, Path] = {}
    for filename, _, _, _ in REFERENCE_FILES:
        # Check canonical name
        p = refs_dir / filename
        if p.exists():
            found[filename] = p
            continue
        # Check alternate names
        for alt in ALTERNATES.get(filename, []):
            p = refs_dir / alt
            if p.exists():
                found[filename] = p
                break
    return found


def setup_reference_files(refs_dir: Path, existing_refs_dir: Path | None = None) -> list[str]:
    """
    Download hg38 reference files or verify existing ones.

    If existing_refs_dir is provided, checks that directory first.
    Only downloads files that are genuinely missing.
    Returns list of failures.
    """
    console.print(Panel("[bold]Step 3 — Reference Files (hg38)[/bold]", style="blue"))
    refs_dir.mkdir(parents=True, exist_ok=True)

    # ── Check existing path first ────────────────────────────────────────────
    search_dir = existing_refs_dir if existing_refs_dir else refs_dir
    if existing_refs_dir:
        console.print(f"  Checking existing references in: [cyan]{existing_refs_dir}[/cyan]")

    existing = _find_existing_refs(search_dir)

    table = Table(show_header=True, header_style="bold cyan")
    table.add_column("File", style="white")
    table.add_column("Size", justify="right", style="yellow")
    table.add_column("Status", justify="center")
    table.add_column("Path", style="dim")

    missing = []
    for filename, gcs_path, description, size_mb in REFERENCE_FILES:
        size_str = f"~{size_mb:,} MB" if size_mb > 1 else "< 1 MB"
        if filename in existing:
            table.add_row(filename, size_str,
                          "[green]✔ found[/green]",
                          str(existing[filename]))
        else:
            table.add_row(filename, size_str,
                          "[red]missing[/red]", "—")
            dest = refs_dir / filename
            missing.append((filename, gcs_path, description, size_mb, dest))

    console.print(table)

    if not missing:
        console.print("  [green]✔[/green]  All reference files already present — skipping download.")
        return []

    total_mb = sum(m[3] for m in missing)
    console.print(
        f"\n  [yellow]{len(missing)} file(s) missing.[/yellow] "
        f"Total download size: ~{total_mb:,} MB"
    )

    if not _ask(f"Download {len(missing)} missing reference file(s) to {refs_dir}?"):
        return [f"Reference file not downloaded: {m[0]}" for m in missing]

    # ── Choose download method ───────────────────────────────────────────────
    use_gsutil = _gsutil_available()
    if not use_gsutil:
        console.print(
            "  [yellow]⚠[/yellow]  gsutil not found. Using wget instead.\n"
            "  [dim]For faster downloads:[/dim] "
            "[cyan]conda install -c conda-forge google-cloud-sdk[/cyan]"
        )

    failures = []
    for filename, gcs_path, description, size_mb, dest in missing:
        console.print(f"\n  Downloading: [bold]{filename}[/bold] ({description})")
        if use_gsutil:
            ok = _download_gsutil(gcs_path, dest)
        else:
            https_url = gcs_path.replace(
                "gs://genomics-public-data",
                "https://storage.googleapis.com/genomics-public-data",
            ).replace(
                "gs://gatk-best-practices",
                "https://storage.googleapis.com/gatk-best-practices",
            )
            ok = _download_wget(https_url, dest)

        if ok:
            console.print(f"  [green]✔[/green]  {filename} downloaded.")
        else:
            console.print(f"  [red]✘[/red]  {filename} download failed.")
            failures.append(f"Download failed: {filename}")

    # ── BWA index ────────────────────────────────────────────────────────────
    fasta = refs_dir / "Homo_sapiens_assembly38.fasta"
    if not fasta.exists() and existing_refs_dir:
        fasta = existing.get("Homo_sapiens_assembly38.fasta", fasta)
    bwa_index_file = fasta.parent / (fasta.name + ".bwt")
    if fasta.exists() and not bwa_index_file.exists():
        console.print("\n  [yellow]⚠[/yellow]  BWA index not found for reference.")
        if _ask("Build BWA index now? (takes ~60 min, ~8 GB disk space)"):
            if not _bwa_index(fasta):
                failures.append("BWA index failed")
            else:
                console.print("  [green]✔[/green]  BWA index built.")

    return failures


def setup_annovar_databases(annovar_bin: Path, annovar_db: Path) -> list[str]:
    """Download missing ANNOVAR databases. Returns list of failures."""
    console.print(Panel("[bold]Step 4 — ANNOVAR Databases[/bold]", style="blue"))

    annotate_pl = annovar_bin / "annotate_variation.pl"
    if not annotate_pl.exists():
        msg = (
            f"ANNOVAR not found at {annotate_pl}\n"
            "  Register and download from: https://annovar.openbioinformatics.org\n"
            f"  Then re-run: exomeflow setup --annovar-bin {annovar_bin} "
            f"--annovar-db {annovar_db}"
        )
        console.print(f"  [red]✘[/red]  {msg}")
        return [msg]

    annovar_db.mkdir(parents=True, exist_ok=True)

    # Show status table
    table = Table(show_header=True, header_style="bold cyan")
    table.add_column("Database", style="white")
    table.add_column("Description")
    table.add_column("Size", justify="right", style="yellow")
    table.add_column("Status", justify="center")

    missing = []
    for db_name, description, size_mb in ANNOVAR_DATABASES:
        db_file = annovar_db / f"hg38_{db_name}.txt"
        exists = db_file.exists()
        status = "[green]exists[/green]" if exists else "[red]missing[/red]"
        size_str = f"~{size_mb:,} MB"
        table.add_row(db_name, description, size_str, status)
        if not exists:
            missing.append((db_name, description, size_mb))

    console.print(table)

    if not missing:
        console.print("  [green]✔[/green]  All ANNOVAR databases present.")
        return []

    total_mb = sum(m[2] for m in missing)
    console.print(
        f"\n  [yellow]{len(missing)} database(s) missing.[/yellow] "
        f"Total size: ~{total_mb:,} MB (~{total_mb // 1024} GB)"
    )

    # Ask per-database or all at once
    if _ask(f"Download all {len(missing)} missing ANNOVAR database(s)?"):
        to_download = missing
    else:
        to_download = []
        for db_name, description, size_mb in missing:
            if _ask(f"  Download {db_name} ({description}, ~{size_mb:,} MB)?"):
                to_download.append((db_name, description, size_mb))

    failures = []
    for db_name, description, size_mb in to_download:
        console.print(f"\n  Downloading: [bold]{db_name}[/bold] ({description}) ...")
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
            console.print(f"  [green]✔[/green]  {db_name} downloaded.")
        else:
            console.print(f"  [red]✘[/red]  {db_name} download failed.")
            failures.append(f"ANNOVAR db download failed: {db_name}")

    skipped = [db for db in missing if db not in to_download]
    for db_name, description, _ in skipped:
        failures.append(f"ANNOVAR db skipped: {db_name} — download manually later")

    return failures


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_setup(
    refs_dir: Path,
    annovar_bin: Path,
    annovar_db: Path,
    existing_refs_dir: Path | None = None,
) -> int:
    """
    Run the full setup sequence.
    Returns number of failures (0 = all good).
    """
    console.print(Panel(
        "[bold white]ExomeFlow — Automated Environment Setup[/bold white]\n"
        "[dim]This will install tools, download references and ANNOVAR databases.[/dim]",
        style="bold blue",
        expand=False,
    ))

    all_failures: list[str] = []

    all_failures += setup_python_packages()
    all_failures += setup_system_tools()
    all_failures += setup_reference_files(refs_dir, existing_refs_dir)
    all_failures += setup_annovar_databases(annovar_bin, annovar_db)

    # ── Summary ──────────────────────────────────────────────────────────
    console.print()
    if not all_failures:
        console.print(Panel(
            "[bold green]✔  Setup complete! All dependencies and databases are ready.[/bold green]\n\n"
            "Run your pipeline with:\n"
            f"[cyan]exomeflow run \\\n"
            f"  --input-dir fastq/ \\\n"
            f"  --reference {refs_dir}/Homo_sapiens_assembly38.fasta \\\n"
            f"  --dbsnp {refs_dir}/dbsnp_146.hg38.vcf.gz \\\n"
            f"  --mills {refs_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \\\n"
            f"  --known-indels {refs_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz \\\n"
            f"  --annovar-bin {annovar_bin} \\\n"
            f"  --annovar-db {annovar_db} \\\n"
            f"  --output results/[/cyan]",
            style="green",
        ))
    else:
        console.print(Panel(
            f"[bold yellow]Setup completed with {len(all_failures)} issue(s):[/bold yellow]\n\n"
            + "\n".join(f"  [red]•[/red] {f}" for f in all_failures),
            style="yellow",
        ))

    return len(all_failures)
