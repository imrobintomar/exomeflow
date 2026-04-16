"""
ExomeFlow — Automated environment setup.

Run via:
    exomeflow setup

What happens:
  1. Auto-detects bundled GATK and ANNOVAR from the ExomeFlow folder
  2. Installs missing system tools (bwa, samtools, fastp, perl) via conda
  3. Asks the user for hg38 reference paths — or downloads them
  4. Asks the user for ANNOVAR humandb path — or downloads databases
  5. Saves everything to ~/.exomeflow/config.json for zero-arg `exomeflow run`
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

console = Console()

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

GATK_VERSION = "4.6.2.0"
GATK_URL = (
    f"https://github.com/broadinstitute/gatk/releases/download/"
    f"{GATK_VERSION}/gatk-{GATK_VERSION}.zip"
)

# Broad Institute GATK resource bundle — publicly accessible (no auth needed)
_GCS_BASE   = "gs://gcp-public-data--broad-references/hg38/v0"
_HTTPS_BASE = "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"

REFERENCE_FILES: list[tuple[str, str, int]] = [
    # (filename, description, approx_size_MB)
    # Core FASTA + pre-built indexes (skip `bwa index` — already in bucket)
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
]

# Alternate filenames users may already have on disk
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
}

ANNOVAR_DATABASES: list[tuple[str, str, int]] = [
    # (db_name, description, approx_size_MB)
    ("refGene",            "Gene annotation (RefSeq)",               80),
    ("avsnp150",           "dbSNP 150 with allele frequencies",  10_000),
    ("clinvar_20240416",   "ClinVar clinical significance",          120),
    ("gnomad41_exome",     "gnomAD v4.1 exome allele frequencies", 8_000),
    ("gnomad41_genome",    "gnomAD v4.1 genome allele frequencies",30_000),
    ("dbnsfp47a",          "dbNSFP 4.7a functional predictions",  50_000),
    ("cosmic84_coding",    "COSMIC v84 somatic mutations",           500),
    ("exac03",             "ExAC 0.3 allele frequencies",          1_200),
]

_SYSTEM_TOOLS: list[tuple[str, str]] = [
    # (binary_name, conda_spec)
    ("bwa",      "bioconda::bwa"),
    ("samtools", "bioconda::samtools"),
    ("fastp",    "bioconda::fastp"),
    ("perl",     "conda-forge::perl"),
]

CONFIG_PATH = Path.home() / ".exomeflow" / "config.json"

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
    Checks ExomeFlow source root first, then PATH and common locations.
    """
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
    ]:
        if candidate.is_file():
            return candidate

    return None


def detect_annovar_bin() -> Path | None:
    """
    Return path to the ANNOVAR directory containing table_annovar.pl.
    Checks ExomeFlow source root first, then common locations.
    """
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


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def _run_silent(cmd: list[str]) -> bool:
    result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return result.returncode == 0


def _run_visible(cmd: list[str]) -> bool:
    result = subprocess.run(cmd, text=True)
    return result.returncode == 0


def _ask(prompt: str, default_yes: bool = False) -> bool:
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
    """Download with wget (resume-capable) or curl."""
    for tool, cmd in [
        ("wget",  ["wget", "-c", "-q", "--show-progress", url, "-O", str(dest)]),
        ("curl",  ["curl", "-L", "--progress-bar", "-C", "-", url, "-o", str(dest)]),
    ]:
        if shutil.which(tool):
            return _run_visible(cmd)
    console.print("  [red]✘[/red]  Neither wget nor curl found.")
    return False


def _gsutil_cp(gcs_path: str, dest: Path) -> bool:
    result = subprocess.run(["gsutil", "-m", "cp", gcs_path, str(dest)], text=True)
    return result.returncode == 0


# ---------------------------------------------------------------------------
# Step 1 — Bundled tools (GATK + ANNOVAR)
# ---------------------------------------------------------------------------

def _step_bundled_tools() -> tuple[Path | None, Path | None]:
    """
    Detect bundled GATK and ANNOVAR. Add GATK to PATH.
    Returns (gatk_path, annovar_bin).
    """
    console.print(Panel("[bold]Step 1 — Bundled Tools[/bold]", style="blue"))

    # GATK
    gatk = detect_gatk_bin()
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
            f"  [red]✘[/red]  GATK not found.\n"
            f"  Place gatk-{GATK_VERSION}/ inside the ExomeFlow folder or download from:\n"
            f"  [cyan]{GATK_URL}[/cyan]"
        )

    # ANNOVAR
    annovar = detect_annovar_bin()
    if annovar:
        console.print(f"  [green]✔[/green]  ANNOVAR: [cyan]{annovar}[/cyan]")
    else:
        console.print(
            "  [red]✘[/red]  ANNOVAR not found.\n"
            "  Place the ANNOVAR folder (with table_annovar.pl) inside the ExomeFlow folder as 'annovar/'.\n"
            "  Register + download: [cyan]https://annovar.openbioinformatics.org[/cyan]"
        )

    return gatk, annovar


# ---------------------------------------------------------------------------
# Step 2 — System tools
# ---------------------------------------------------------------------------

def _step_system_tools() -> list[str]:
    """Install missing system tools via conda. Returns list of failures."""
    console.print(Panel("[bold]Step 2 — System Tools[/bold]", style="blue"))
    failures = []

    conda = shutil.which("conda") or shutil.which("mamba")
    if not conda:
        console.print(
            "  [yellow]⚠[/yellow]  conda/mamba not found — cannot auto-install tools.\n"
            "  Install Miniconda: [cyan]https://docs.conda.io/en/latest/miniconda.html[/cyan]"
        )
        return ["conda not found — install tools manually"]

    for binary, spec in _SYSTEM_TOOLS:
        if shutil.which(binary):
            console.print(f"  [green]✔[/green]  {binary} already installed")
            continue
        console.print(f"  [yellow]→[/yellow]  Installing {binary} ...")
        channel, pkg = spec.split("::")
        ok = _run_silent([conda, "install", "-y", "-c", channel, pkg])
        if ok:
            console.print(f"  [green]✔[/green]  {binary} installed")
        else:
            console.print(f"  [red]✘[/red]  {binary} install failed — try: conda install -c {channel} {pkg}")
            failures.append(f"{binary} install failed")

    return failures


# ---------------------------------------------------------------------------
# Step 3 — Reference files
# ---------------------------------------------------------------------------

def _scan_refs(directory: Path) -> dict[str, Path]:
    """Return {canonical_name: found_path} for all reference files found in directory."""
    found: dict[str, Path] = {}
    for filename, _, _ in REFERENCE_FILES:
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


def _step_reference_files(refs_dir: Path, existing_refs_dir: Path | None) -> tuple[Path | None, list[str]]:
    """
    Locate or download hg38 reference files.
    Returns (resolved_refs_dir, failures).
    """
    console.print(Panel("[bold]Step 3 — Reference Genome (hg38)[/bold]", style="blue"))
    failures: list[str] = []

    # If user passed --existing-refs, check there first
    if existing_refs_dir:
        found = _scan_refs(existing_refs_dir)
        if found:
            console.print(f"  [green]✔[/green]  Found {len(found)}/{len(REFERENCE_FILES)} reference files in [cyan]{existing_refs_dir}[/cyan]")
            return existing_refs_dir, []
        console.print(f"  [yellow]⚠[/yellow]  No reference files found in {existing_refs_dir}")

    # Check default location
    found = _scan_refs(refs_dir)
    if len(found) >= 4:  # the 4 main files (fasta + 3 VCFs)
        console.print(f"  [green]✔[/green]  Reference files already present in [cyan]{refs_dir}[/cyan]")
        return refs_dir, []

    # Ask user
    console.print(
        "  hg38 reference files are required (~13 GB total).\n"
        "  These include: reference FASTA, dbSNP, Mills indels, known indels."
    )

    if _ask("Do you already have hg38 reference files on this machine?"):
        user_path = _ask_path("Enter the path to the directory containing your reference files")
        if user_path:
            found = _scan_refs(user_path)
            if found:
                console.print(f"  [green]✔[/green]  Found {len(found)} reference file(s) in [cyan]{user_path}[/cyan]")
                return user_path, []
            else:
                console.print(f"  [red]✘[/red]  No recognised reference files found in {user_path}")
                failures.append(f"No reference files found in {user_path}")
                return None, failures
        else:
            console.print("  [red]✘[/red]  Path not found or not entered.")
            failures.append("Reference path not provided")
            return None, failures

    # Offer download
    total_mb = sum(s for _, _, s in REFERENCE_FILES)
    if not _ask(f"Download reference files to {refs_dir}? (~{total_mb // 1024} GB, may take hours)", default_yes=True):
        failures.append("Reference files not downloaded — required for pipeline")
        return None, failures

    refs_dir.mkdir(parents=True, exist_ok=True)
    use_gsutil = bool(shutil.which("gsutil"))
    if not use_gsutil:
        console.print(
            "  [dim]gsutil not found — using wget. "
            "For faster downloads: conda install -c conda-forge google-cloud-sdk[/dim]"
        )

    for filename, description, size_mb in REFERENCE_FILES:
        dest = refs_dir / filename
        if dest.exists():
            console.print(f"  [green]✔[/green]  {filename} already present")
            continue
        console.print(f"  [cyan]→[/cyan]  Downloading {filename} ({description}, ~{size_mb:,} MB) ...")
        if use_gsutil:
            ok = _gsutil_cp(f"{_GCS_BASE}/{filename}", dest)
        else:
            ok = _download_file(f"{_HTTPS_BASE}/{filename}", dest)
        if ok:
            console.print(f"  [green]✔[/green]  {filename}")
        else:
            console.print(f"  [red]✘[/red]  {filename} download failed")
            failures.append(f"Download failed: {filename}")

    return refs_dir if not failures else None, failures


# ---------------------------------------------------------------------------
# Step 4 — ANNOVAR databases
# ---------------------------------------------------------------------------

def _step_annovar_databases(annovar_bin: Path, default_db: Path) -> tuple[Path | None, list[str]]:
    """
    Locate or download ANNOVAR databases.
    Returns (humandb_path, failures).
    """
    console.print(Panel("[bold]Step 4 — ANNOVAR Annotation Databases[/bold]", style="blue"))
    failures: list[str] = []

    if not annovar_bin:
        console.print("  [yellow]⚠[/yellow]  ANNOVAR not detected — skipping database step.")
        failures.append("ANNOVAR not found — databases not set up")
        return None, failures

    annotate_pl = annovar_bin / "annotate_variation.pl"

    # Check existing db in the annovar_bin/humandb location
    for candidate_db in [default_db, annovar_bin / "humandb"]:
        if candidate_db.exists():
            existing = [d for d, _, _ in ANNOVAR_DATABASES
                        if (candidate_db / f"hg38_{d}.txt").exists()]
            if existing:
                console.print(
                    f"  [green]✔[/green]  Found {len(existing)}/{len(ANNOVAR_DATABASES)} "
                    f"databases in [cyan]{candidate_db}[/cyan]"
                )
                return candidate_db, []

    # Ask user
    total_gb = sum(s for _, _, s in ANNOVAR_DATABASES) // 1024
    console.print(
        f"  ANNOVAR databases are required for variant annotation (~{total_gb} GB total).\n"
        "  These include: refGene, ClinVar, gnomAD, dbNSFP, COSMIC, avSNP."
    )

    if _ask("Do you already have an ANNOVAR humandb directory?"):
        user_path = _ask_path("Enter the path to your humandb directory")
        if user_path and user_path.is_dir():
            existing = [d for d, _, _ in ANNOVAR_DATABASES
                        if (user_path / f"hg38_{d}.txt").exists()]
            console.print(
                f"  [green]✔[/green]  Found {len(existing)}/{len(ANNOVAR_DATABASES)} "
                f"databases in [cyan]{user_path}[/cyan]"
            )
            return user_path, []
        else:
            console.print("  [red]✘[/red]  Path not found or not entered.")
            failures.append("ANNOVAR humandb path not valid")
            return None, failures

    # Offer download
    db_dir = default_db
    if not _ask(f"Download ANNOVAR databases to {db_dir}? (~{total_gb} GB)", default_yes=True):
        failures.append("ANNOVAR databases not downloaded — required for annotation")
        return None, failures

    db_dir.mkdir(parents=True, exist_ok=True)

    for db_name, description, size_mb in ANNOVAR_DATABASES:
        db_file = db_dir / f"hg38_{db_name}.txt"
        if db_file.exists():
            console.print(f"  [green]✔[/green]  {db_name} already present")
            continue
        console.print(f"  [cyan]→[/cyan]  Downloading {db_name} ({description}, ~{size_mb:,} MB) ...")
        result = subprocess.run(
            ["perl", str(annotate_pl), "-buildver", "hg38",
             "-downdb", "-webfrom", "annovar", db_name, str(db_dir)],
            text=True,
        )
        if result.returncode == 0:
            console.print(f"  [green]✔[/green]  {db_name}")
        else:
            console.print(f"  [red]✘[/red]  {db_name} download failed")
            failures.append(f"ANNOVAR db download failed: {db_name}")

    return db_dir if not failures else None, failures


# ---------------------------------------------------------------------------
# Pre-flight dependency check (called on every `exomeflow run`)
# ---------------------------------------------------------------------------

_TOOL_CHECKS: list[tuple[str, list[str], str]] = [
    # (name, version_cmd, min_version_string)
    ("bwa",      ["bwa"],                    "BWA aligner"),
    ("samtools", ["samtools", "--version"],   "SAMtools"),
    ("fastp",    ["fastp", "--version"],      "fastp QC"),
    ("perl",     ["perl", "--version"],       "Perl interpreter"),
    ("gatk",     ["gatk", "--version"],       "GATK 4"),
]

_REF_KEYS: list[tuple[str, str]] = [
    # (config_key, label)
    ("reference",    "Reference FASTA (hg38)"),
    ("dbsnp",        "dbSNP VCF"),
    ("mills",        "Mills indels VCF"),
    ("known_indels", "Known indels VCF"),
]

_ANNOVAR_KEYS: list[tuple[str, str]] = [
    ("annovar_bin", "ANNOVAR scripts directory"),
    ("annovar_db",  "ANNOVAR humandb directory"),
]


def check_and_fix_dependencies() -> dict:
    """
    Run a pre-flight check before every pipeline run.

    Prints a status table showing every dependency.
    If anything is missing: runs the relevant fix step automatically.
    Returns the final config dict (with all paths resolved) on success.
    Raises SystemExit(1) if any required dependency cannot be fixed.
    """
    from rich.panel import Panel

    cfg = load_config()

    # ── Build status for all checks ─────────────────────────────────────────
    issues: list[str] = []   # categories that need fixing

    table = Table(show_header=True, header_style="bold cyan", box=None, padding=(0, 2))
    table.add_column("Category",    style="bold white",  min_width=22)
    table.add_column("Dependency",  style="white",       min_width=30)
    table.add_column("Status",      justify="center",    min_width=12)

    # Tools
    tool_missing: list[str] = []
    for name, _cmd, label in _TOOL_CHECKS:
        ok = bool(shutil.which(name))
        status = "[green]✔ found[/green]" if ok else "[red]✘ missing[/red]"
        table.add_row("Tools" if name == "bwa" else "", f"{label} ({name})", status)
        if not ok:
            tool_missing.append(name)
    if tool_missing:
        issues.append("tools")

    # GATK via bundled path even when not on PATH
    gatk_path = detect_gatk_bin()
    if gatk_path and not shutil.which("gatk"):
        # Fix PATH inline — update the "gatk" row retroactively isn't easy,
        # but we'll show it as OK in the table since we can resolve it
        tool_missing = [t for t in tool_missing if t != "gatk"]
        if not tool_missing:
            issues = [i for i in issues if i != "tools"]

    # Reference files
    ref_missing: list[str] = []
    for key, label in _REF_KEYS:
        val = cfg.get(key)
        ok = bool(val and Path(val).exists())
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
        ok = bool(val and Path(val).exists())
        status = "[green]✔ found[/green]" if ok else "[red]✘ missing[/red]"
        table.add_row("ANNOVAR" if key == "annovar_bin" else "", label, status)
        if not ok:
            annovar_missing.append(key)
    if annovar_missing:
        issues.append("annovar")

    # ── Print the table ──────────────────────────────────────────────────────
    console.print()
    console.print(Panel(table, title="[bold]ExomeFlow — Dependency Check[/bold]",
                        border_style="blue", expand=False))

    if not issues:
        console.print("  [green]✔  All dependencies satisfied.[/green] Starting pipeline...\n")
        return cfg

    # ── Fix what's missing ──────────────────────────────────────────────────
    console.print(f"\n  [yellow]⚠  Fixing {len(issues)} missing requirement(s)...[/yellow]\n")

    new_cfg: dict = {}

    # Fix tools
    if "tools" in issues:
        console.print(Panel("[bold]Fixing: System Tools[/bold]", style="yellow"))
        # GATK — add bundled path to PATH
        if "gatk" in tool_missing:
            if gatk_path:
                _add_to_path(gatk_path)
                console.print(f"  [green]✔[/green]  GATK resolved: [cyan]{gatk_path}[/cyan]")
                tool_missing.remove("gatk")
                new_cfg["gatk_bin"] = gatk_path
            else:
                console.print(
                    f"  [red]✘[/red]  GATK not found. Place [cyan]gatk-{GATK_VERSION}/[/cyan] "
                    "inside the ExomeFlow folder."
                )
        # Other tools via conda
        remaining_tools = [t for t in tool_missing if t != "gatk"]
        if remaining_tools:
            _step_system_tools()

    # Fix reference files
    if "refs" in issues:
        console.print(Panel("[bold]Fixing: Reference Files[/bold]", style="yellow"))
        refs_dir = Path(cfg.get("refs_dir", Path.home() / ".exomeflow" / "refs"))
        resolved_refs, ref_failures = _step_reference_files(refs_dir, existing_refs_dir=None)
        if resolved_refs is None:
            console.print("  [red]✘[/red]  Reference files are required. Cannot continue.")
            raise SystemExit(1)
        new_cfg["refs_dir"]     = resolved_refs
        new_cfg["reference"]    = _find_ref(resolved_refs, "Homo_sapiens_assembly38.fasta")
        new_cfg["dbsnp"]        = _find_ref(resolved_refs, "Homo_sapiens_assembly38.dbsnp138.vcf.gz")
        new_cfg["mills"]        = _find_ref(resolved_refs, "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
        new_cfg["known_indels"] = _find_ref(resolved_refs, "Homo_sapiens_assembly38.known_indels.vcf.gz")

    # Fix ANNOVAR
    if "annovar" in issues:
        console.print(Panel("[bold]Fixing: ANNOVAR[/bold]", style="yellow"))
        annovar_bin = detect_annovar_bin()
        if annovar_bin is None:
            console.print(
                "  [red]✘[/red]  ANNOVAR not found. Place the [cyan]annovar/[/cyan] folder "
                "inside the ExomeFlow directory.\n"
                "  Register + download: https://annovar.openbioinformatics.org"
            )
            raise SystemExit(1)
        default_db = Path(cfg.get("annovar_db", annovar_bin / "humandb"))
        resolved_db, db_failures = _step_annovar_databases(annovar_bin, default_db)
        if resolved_db is None:
            console.print("  [red]✘[/red]  ANNOVAR databases are required. Cannot continue.")
            raise SystemExit(1)
        new_cfg["annovar_bin"] = annovar_bin
        new_cfg["annovar_db"]  = resolved_db

    # Save any newly resolved paths
    if new_cfg:
        save_config(new_cfg)
        cfg = load_config()
        console.print(f"\n  [dim]Configuration updated → {CONFIG_PATH}[/dim]")

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
# Main entry point
# ---------------------------------------------------------------------------

def run_setup(
    refs_dir: Path,
    annovar_bin: Path | None = None,
    annovar_db: Path | None = None,
    existing_refs_dir: Path | None = None,
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

    # ── Step 1: Bundled GATK + ANNOVAR ──────────────────────────────────────
    gatk_path, resolved_annovar_bin = _step_bundled_tools()
    if gatk_path:
        config["gatk_bin"] = gatk_path
    effective_annovar_bin = annovar_bin or resolved_annovar_bin

    # ── Step 2: System tools ────────────────────────────────────────────────
    all_failures += _step_system_tools()

    # ── Step 3: Reference files ─────────────────────────────────────────────
    resolved_refs, ref_failures = _step_reference_files(refs_dir, existing_refs_dir)
    all_failures += ref_failures
    if resolved_refs:
        # Resolve actual paths — handle both dbsnp138 and alternate names
        def _find_in(d: Path, canonical: str) -> Path:
            for alt in [canonical] + _REF_ALTERNATES.get(canonical, []):
                p = d / alt
                if p.exists():
                    return p
            return d / canonical  # fallback even if not present yet

        config["refs_dir"]     = resolved_refs
        config["reference"]    = _find_in(resolved_refs, "Homo_sapiens_assembly38.fasta")
        config["dbsnp"]        = _find_in(resolved_refs, "Homo_sapiens_assembly38.dbsnp138.vcf.gz")
        config["mills"]        = _find_in(resolved_refs, "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
        config["known_indels"] = _find_in(resolved_refs, "Homo_sapiens_assembly38.known_indels.vcf.gz")

    # ── Step 4: ANNOVAR databases ────────────────────────────────────────────
    if effective_annovar_bin:
        default_db = annovar_db or (effective_annovar_bin / "humandb")
        resolved_db, db_failures = _step_annovar_databases(effective_annovar_bin, default_db)
        all_failures += db_failures
        if resolved_db:
            config["annovar_bin"] = effective_annovar_bin
            config["annovar_db"]  = resolved_db

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
