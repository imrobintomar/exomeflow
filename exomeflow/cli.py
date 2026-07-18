"""
Typer-based CLI for ExomeFlow.

Entry point registered in pyproject.toml:
    exomeflow = "exomeflow.cli:app"
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from exomeflow import __author__, __email__, __version__

_HELP = """
[bold cyan]ExomeFlow[/bold cyan] — Production Whole Exome/Genome Sequencing pipeline.

[bold]Quick start[/bold]

  [dim]# Install[/dim]
  pip install exomeflow

  [dim]# Run — first launch walks you through one-time setup automatically[/dim]
  exomeflow run --input-dir fastq/ --output results/

[bold]Commands[/bold]
  [green]run[/green]     Execute the full WES analysis pipeline on paired FASTQ files
  [green]setup[/green]   Re-run setup (change reference paths, download new databases)
  [green]doctor[/green]  Read-only pre-flight report of what's found/missing, before you run

[bold]What happens on first run[/bold]
  1. GATK is auto-detected or auto-downloaded. ANNOVAR is auto-detected, or
     you're prompted for where you extracted it (it requires free personal
     registration before download — the one thing that can't be automated;
     run [green]exomeflow doctor[/green] any time to see exactly what's missing)
  2. Missing tools (bwa, samtools, fastp, perl) are installed via conda/mamba
     — or a self-bootstrapped micromamba if neither is already installed
  3. You are asked for hg38/GRCh37 reference paths (or they are downloaded)
  4. You are asked for ANNOVAR database paths (or they are downloaded)
  5. HPO, InterVar (ACMG), and MultiQC are provisioned automatically too
  6. Threads/JVM memory are sized from your CPU/RAM automatically
  7. Everything is saved — future runs need no extra arguments
  (add [green]--yes[/green] to skip all prompts for unattended/CI runs)

[bold]Beyond default germline single-sample calling[/bold]
  --joint-genotyping   Cohort mode: one shared VCF instead of per-sample files
  --mode somatic       Tumor-only Mutect2 instead of HaplotypeCaller
  --cnv                Read-depth CNV calling alongside SNP/INDEL calling
  --genome-build       hg38 (default) or GRCh37
"""

app = typer.Typer(
    name="exomeflow",
    help=_HELP,
    add_completion=False,
    rich_markup_mode="rich",
    no_args_is_help=True,
)

console = Console()


def _version_callback(value: bool) -> None:
    if value:
        console.print(f"ExomeFlow v{__version__}")
        console.print(f"{__author__}, AIIMS New Delhi <{__email__}>")
        raise typer.Exit()


_VALID_GENOME_BUILDS = ("hg38", "GRCh37")


def _validate_genome_build(genome_build: Optional[str]) -> None:
    """Shared by run/setup — was duplicated verbatim in both commands."""
    if genome_build is not None and genome_build not in _VALID_GENOME_BUILDS:
        console.print(
            f"  [red]✘[/red]  --genome-build must be 'hg38' or 'GRCh37', got '{genome_build}'."
        )
        raise typer.Exit(code=1)


@app.callback()
def _main(
    version: Optional[bool] = typer.Option(
        None,
        "--version", "-v",
        callback=_version_callback,
        is_eager=True,
        help="Print version and exit.",
    ),
) -> None:
    pass


def _ensure_ready(genome_build: Optional[str], call_cnv: bool, assume_yes: bool, mode: str = "germline") -> dict:
    """
    Called at the start of every `exomeflow run`.
    Checks all dependencies, fixes what's missing, returns the resolved config.
    """
    from exomeflow.setup_env import check_and_fix_dependencies
    try:
        return check_and_fix_dependencies(
            genome_build=genome_build, call_cnv=call_cnv, assume_yes=assume_yes, mode=mode
        )
    except SystemExit as exc:
        raise typer.Exit(code=int(exc.code or 1))


@app.command("run")
def run_command(
    input_dir: Path = typer.Option(
        ...,
        "--input-dir", "-i",
        help="Directory containing paired FASTQ files (*_1.fastq.gz/*_2.fastq.gz or *_R1.fastq.gz/*_R2.fastq.gz).",
        exists=True, file_okay=False, dir_okay=True, readable=True,
    ),
    output: Path = typer.Option(
        Path("results"),
        "--output", "-o",
        help="Root output directory (created if absent).",
    ),
    # ── Reference overrides (optional — auto-resolved from saved config) ──
    reference: Optional[Path] = typer.Option(None, "--reference", "-r",
        help="Reference genome FASTA (hg38 or GRCh37 — see --genome-build). "
             "Auto-resolved from saved config if omitted."),
    dbsnp: Optional[Path] = typer.Option(None, "--dbsnp",
        help="dbSNP VCF. Auto-resolved from saved config if omitted."),
    mills: Optional[Path] = typer.Option(None, "--mills",
        help="Mills indels VCF. Auto-resolved from saved config if omitted."),
    known_indels: Optional[Path] = typer.Option(None, "--known-indels",
        help="Known indels VCF. Auto-resolved from saved config if omitted."),
    annovar_bin: Optional[Path] = typer.Option(None, "--annovar-bin",
        help="ANNOVAR directory. Auto-resolved from saved config if omitted."),
    annovar_db: Optional[Path] = typer.Option(None, "--annovar-db",
        help="ANNOVAR humandb directory. Auto-resolved from saved config if omitted."),
    annovar_protocols: Optional[str] = typer.Option(None, "--annovar-protocols",
        help="Override the ANNOVAR --protocol list (comma-separated db names). "
             "Use this if your existing humandb has different database versions "
             "than the shipped default — e.g. a newer/older ClinVar build."),
    annovar_operations: Optional[str] = typer.Option(None, "--annovar-operations",
        help="Override the ANNOVAR --operation list. Must have the same number "
             "of comma-separated entries as --annovar-protocols."),
    # ── Pipeline options ──────────────────────────────────────────────────
    intervals: Optional[Path] = typer.Option(None, "--intervals",
        help="Exome capture BED file. Omit for whole-genome mode."),
    interval_padding: int = typer.Option(100, "--interval-padding",
        help="Base-pair padding around each target interval."),
    threads: Optional[int] = typer.Option(None, "--threads", "-t",
        help="Threads for BWA MEM and GATK HaplotypeCaller. Auto-detected from "
             "CPU count if omitted.", min=1),
    fastp_threads: int = typer.Option(8, "--fastp-threads",
        help="Threads for fastp.", min=1),
    annovar_threads: int = typer.Option(24, "--annovar-threads",
        help="Threads for ANNOVAR.", min=1),
    max_workers: int = typer.Option(1, "--max-workers",
        help="Number of samples to process in parallel.", min=1),
    java_opts: Optional[str] = typer.Option(None, "--java-opts",
        help="JVM options passed via JAVA_OPTS. Auto-sized from available RAM "
             "if omitted (60% of available, 4-80g)."),
    # ── V2 mode flags (all default to v1 behavior) ─────────────────────────
    mode: str = typer.Option("germline", "--mode",
        help="Variant-calling mode: 'germline' (HaplotypeCaller) or 'somatic' "
             "(tumor-only Mutect2)."),
    genome_build: Optional[str] = typer.Option(None, "--genome-build",
        help="Reference genome build: 'hg38' or 'GRCh37'. If omitted: uses the "
             "saved choice, or asks interactively, or defaults to hg38 under --yes."),
    joint_genotyping: bool = typer.Option(False, "--joint-genotyping",
        help="Cohort mode: joint-genotype all samples into one shared VCF "
             "instead of per-sample annotated files. Opt-in only."),
    call_cnv: bool = typer.Option(False, "--cnv",
        help="Also call read-depth CNVs per sample (requires --intervals)."),
    germline_resource: Optional[Path] = typer.Option(None, "--germline-resource",
        help="gnomAD AF-only VCF for Mutect2 (--mode somatic). Auto-downloaded "
             "(GATK's public resource) if omitted — pass this to use your own instead."),
    panel_of_normals: Optional[Path] = typer.Option(None, "--panel-of-normals",
        help="Panel of Normals VCF for Mutect2 (--mode somatic). Auto-downloaded "
             "(GATK's public 1000 Genomes PoN) if omitted — pass this to use your own "
             "instead. Filters recurrent sequencing artifacts a population AF "
             "resource alone won't catch."),
    assume_yes: bool = typer.Option(False, "--yes", "-y",
        help="Non-interactive: auto-confirm every setup prompt (downloads, etc.) "
             "instead of asking. Needed for unattended/background/CI runs."),
) -> None:
    """
    Run the complete WES/WGS analysis pipeline.

    \b
    On first run, ExomeFlow sets itself up automatically:
      • Detects or auto-downloads GATK and ANNOVAR
      • Installs missing tools (bwa, samtools, fastp, perl) via conda
      • Asks for reference genome paths (or downloads them)
      • Asks for ANNOVAR database paths (or downloads them)
      • Auto-provisions HPO, InterVar (ACMG), and MultiQC
      • Sizes --threads/--java-opts from your CPU/RAM automatically
      • Saves everything — future runs need no extra arguments

    \b
    Default workflow (germline, single-sample)
    -------------------------------------------
    FASTQ → fastp → BWA MEM → BAM processing → BQSR → HaplotypeCaller
         → hard-filter → ANNOVAR → HPO terms → ACMG classification
    ...then once, across all samples: MultiQC rollup.
    Every sample gets its own separate annotated output file by default.

    \b
    Examples
    --------
    exomeflow run --input-dir fastq/ --output results/
    exomeflow run --input-dir fastq/ --output results/ --mode somatic
    exomeflow run --input-dir fastq/ --output results/ --joint-genotyping --intervals targets.bed
    """
    from exomeflow.config import Config, intervals_present
    from exomeflow.pipeline import run_pipeline
    from exomeflow.setup_env import detect_gatk_bin

    if mode not in ("germline", "somatic"):
        console.print(f"  [red]✘[/red]  --mode must be 'germline' or 'somatic', got '{mode}'.")
        raise typer.Exit(code=1)
    _validate_genome_build(genome_build)
    if (joint_genotyping or call_cnv) and not intervals_present(intervals):
        console.print(
            "  [red]✘[/red]  --joint-genotyping and --cnv both require --intervals "
            "(a bounded region is needed before either can run)."
        )
        raise typer.Exit(code=1)
    if (annovar_protocols is None) != (annovar_operations is None):
        console.print(
            "  [red]✘[/red]  --annovar-protocols and --annovar-operations must be "
            "given together (ANNOVAR pairs each protocol with an operation)."
        )
        raise typer.Exit(code=1)
    if annovar_protocols is not None and (
        annovar_protocols.count(",") != annovar_operations.count(",")
    ):
        console.print(
            "  [red]✘[/red]  --annovar-protocols and --annovar-operations have a "
            "different number of comma-separated entries."
        )
        raise typer.Exit(code=1)

    # ── Fail fast on a bad --input-dir before the setup wizard (which can
    # take hours on first run) ever starts — nothing is more frustrating
    # than a multi-hour download finishing only to hit "no FASTQ files". ──
    from exomeflow.utils import detect_samples
    try:
        detect_samples(input_dir)
    except FileNotFoundError as exc:
        console.print(f"  [red]✘[/red]  {exc}")
        raise typer.Exit(code=1)

    # ── Auto-size --threads/--java-opts from actual system resources
    # instead of a static default that can either undersubscribe a big box
    # or ask for more RAM than a small one has. ───────────────────────────
    if threads is None or java_opts is None:
        from exomeflow.utils import detect_system_resources, recommend_java_opts, recommend_threads
        resources = detect_system_resources(output)
        if threads is None:
            threads = recommend_threads(resources)
        if java_opts is None:
            java_opts = recommend_java_opts(resources)
        console.print(
            f"  [dim]Detected {resources.cpu_count} CPUs, "
            f"{resources.available_ram_gb:.0f}/{resources.total_ram_gb:.0f} GB RAM available, "
            f"{resources.free_disk_gb:.0f} GB free disk — "
            f"using --threads {threads} --java-opts \"{java_opts}\"[/dim]"
        )

    # ── First-run setup (skipped if config already complete) ─────────────
    # genome_build resolution: explicit flag > saved config > interactive
    # prompt > hg38 default under --yes (see setup_env._ask_genome_build).
    # check_and_fix_dependencies() unconditionally sets cfg["genome_build"]
    # on every path, so it's always present here — no fallback needed.
    saved = _ensure_ready(genome_build=genome_build, call_cnv=call_cnv, assume_yes=assume_yes, mode=mode)
    genome_build = saved["genome_build"]

    # ── Add GATK to PATH if not already there ────────────────────────────
    if not shutil.which("gatk"):
        gatk_path = Path(saved["gatk_bin"]) if "gatk_bin" in saved else detect_gatk_bin()
        if gatk_path and gatk_path.is_file():
            gatk_dir = str(gatk_path.parent)
            os.environ["PATH"] = gatk_dir + os.pathsep + os.environ.get("PATH", "")
            try:
                gatk_path.chmod(gatk_path.stat().st_mode | 0o111)
            except Exception:
                pass

    # ── Resolve paths: explicit CLI arg → saved config ────────────────────
    def _r(provided: Optional[Path], key: str) -> Path:
        if provided is not None:
            return provided
        return Path(saved[key])

    # Same idea, but optional: germline_resource/panel_of_normals are only
    # ever in `saved` for --mode somatic runs (auto-downloaded by
    # check_and_fix_dependencies), so a plain _r() would KeyError on any
    # germline run — fall back to None instead of requiring the key exist.
    def _r_optional(provided: Optional[Path], key: str) -> Optional[Path]:
        if provided is not None:
            return provided
        val = saved.get(key)
        return Path(val) if val else None

    cfg_overrides: dict = {}
    if annovar_protocols is not None:
        cfg_overrides["annovar_protocols"] = annovar_protocols
    if annovar_operations is not None:
        cfg_overrides["annovar_operations"] = annovar_operations

    cfg = Config(
        input_dir=input_dir,
        output_dir=output,
        reference=_r(reference, "reference"),
        dbsnp=_r(dbsnp, "dbsnp"),
        mills=_r(mills, "mills"),
        known_indels=_r(known_indels, "known_indels"),
        intervals=intervals,
        interval_padding=interval_padding,
        annovar_bin=_r(annovar_bin, "annovar_bin"),
        annovar_db=_r(annovar_db, "annovar_db"),
        threads=threads,
        fastp_threads=fastp_threads,
        annovar_threads=annovar_threads,
        max_workers=max_workers,
        java_opts=java_opts,
        mode=mode,
        genome_build=genome_build,
        joint_genotyping=joint_genotyping,
        call_cnv=call_cnv,
        germline_resource=_r_optional(germline_resource, "germline_resource"),
        panel_of_normals=_r_optional(panel_of_normals, "panel_of_normals"),
        **cfg_overrides,
    )

    failed = run_pipeline(cfg)
    raise typer.Exit(code=min(failed, 1))


@app.command("setup")
def setup_command(
    refs_dir: Path = typer.Option(
        Path.home() / ".exomeflow" / "refs",
        "--refs-dir",
        help="Directory for reference genome files. Default: ~/.exomeflow/refs",
    ),
    existing_refs: Optional[Path] = typer.Option(
        None, "--existing-refs",
        help="Path to existing reference files — skips download if found here.",
    ),
    annovar_bin: Optional[Path] = typer.Option(
        None, "--annovar-bin",
        help="ANNOVAR directory. Auto-detected (saved config, then common "
             "locations) if omitted — you'll be prompted for it if not found.",
    ),
    annovar_db: Optional[Path] = typer.Option(
        None, "--annovar-db",
        help="ANNOVAR humandb directory. Defaults to <annovar-bin>/humandb.",
    ),
    genome_build: Optional[str] = typer.Option(
        None, "--genome-build",
        help="Reference genome build to fetch: 'hg38' or 'GRCh37'. If omitted: "
             "uses the saved choice, or asks interactively, or defaults to hg38 under --yes.",
    ),
    assume_yes: bool = typer.Option(False, "--yes", "-y",
        help="Non-interactive: auto-confirm every setup prompt instead of asking."),
) -> None:
    """
    Re-run setup: change reference paths, download new databases, or repair config.

    \b
    You do NOT need to run this before your first `exomeflow run`.
    First-time setup happens automatically when you run the pipeline.

    \b
    Use this command to:
      • Switch to different reference files
      • Download additional ANNOVAR databases
      • Reset and repair the saved configuration

    \b
    Examples
    --------
    # Re-run setup (uses saved/auto-detected paths):
    exomeflow setup

    # Point to existing reference files (skip download):
    exomeflow setup --existing-refs /data/hg38

    # Fully explicit:
    exomeflow setup \\
      --refs-dir    /data/refs \\
      --annovar-bin /opt/annovar \\
      --annovar-db  /opt/annovar/humandb
    """
    from exomeflow.setup_env import run_setup

    _validate_genome_build(genome_build)

    # ANNOVAR resolution (auto-detect, then an interactive path prompt if not
    # found) happens inside run_setup() itself now, so the rest of setup
    # (GATK, system tools, reference files) still completes even when
    # ANNOVAR isn't resolved yet — it doesn't have to block everything else.
    failed = run_setup(
        refs_dir=refs_dir,
        annovar_bin=annovar_bin,
        annovar_db=annovar_db,
        existing_refs_dir=existing_refs,
        genome_build=genome_build,
        assume_yes=assume_yes,
    )
    raise typer.Exit(code=min(failed, 1))


@app.command("doctor")
def doctor_command(
    genome_build: Optional[str] = typer.Option(
        None, "--genome-build",
        help="Check ANNOVAR database completeness against this build. "
             "Defaults to the saved choice, or hg38.",
    ),
) -> None:
    """
    Pre-flight report: what's found, what's missing, and what will/won't
    resolve itself automatically.

    \b
    Read-only — makes no network requests, installs nothing, and never
    writes ~/.exomeflow/config.json. Safe to run any time, including
    before your first `exomeflow run`, to see the whole picture upfront
    instead of hitting each gap one at a time mid-setup.

    \b
    Examples
    --------
    exomeflow doctor
    exomeflow doctor --genome-build GRCh37
    """
    from exomeflow.setup_env import run_doctor

    _validate_genome_build(genome_build)
    run_doctor(genome_build=genome_build)
