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

from exomeflow import __version__

_HELP = """
[bold cyan]ExomeFlow[/bold cyan] — Production Whole Exome Sequencing pipeline.

[bold]Quick start[/bold]

  [dim]# Install[/dim]
  pip install exomeflow

  [dim]# Run — first launch walks you through one-time setup automatically[/dim]
  exomeflow run --input-dir fastq/ --output results/

[bold]Commands[/bold]
  [green]run[/green]     Execute the full WES analysis pipeline on paired FASTQ files
  [green]setup[/green]   Re-run setup (change reference paths, download new databases)

[bold]What happens on first run[/bold]
  1. Bundled GATK and ANNOVAR are detected automatically
  2. Missing tools (bwa, samtools, fastp, perl) are installed via conda
  3. You are asked for reference genome paths (or they are downloaded)
  4. You are asked for ANNOVAR database paths (or they are downloaded)
  5. Everything is saved — future runs need no extra arguments
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
        raise typer.Exit()


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


def _ensure_ready() -> dict:
    """
    Called at the start of every `exomeflow run`.
    Checks all dependencies, fixes what's missing, returns the resolved config.
    """
    from exomeflow.setup_env import check_and_fix_dependencies
    try:
        return check_and_fix_dependencies()
    except SystemExit as exc:
        raise typer.Exit(code=int(exc.code or 1))


@app.command("run")
def run_command(
    input_dir: Path = typer.Option(
        ...,
        "--input-dir", "-i",
        help="Directory containing paired FASTQ files (*_1.fastq.gz / *_2.fastq.gz).",
        exists=True, file_okay=False, dir_okay=True, readable=True,
    ),
    output: Path = typer.Option(
        Path("results"),
        "--output", "-o",
        help="Root output directory (created if absent).",
    ),
    # ── Reference overrides (optional — auto-resolved from saved config) ──
    reference: Optional[Path] = typer.Option(None, "--reference", "-r",
        help="hg38 FASTA. Auto-resolved from saved config if omitted."),
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
    # ── Pipeline options ──────────────────────────────────────────────────
    intervals: Path = typer.Option(Path(""), "--intervals",
        help="Exome capture BED file. Omit for whole-genome mode."),
    interval_padding: int = typer.Option(100, "--interval-padding",
        help="Base-pair padding around each target interval."),
    threads: int = typer.Option(24, "--threads", "-t",
        help="Threads for BWA MEM and GATK HaplotypeCaller.", min=1),
    fastp_threads: int = typer.Option(8, "--fastp-threads",
        help="Threads for fastp.", min=1),
    annovar_threads: int = typer.Option(24, "--annovar-threads",
        help="Threads for ANNOVAR.", min=1),
    max_workers: int = typer.Option(1, "--max-workers",
        help="Number of samples to process in parallel.", min=1),
    java_opts: str = typer.Option("-Xmx80g", "--java-opts",
        help="JVM options passed via JAVA_OPTS environment variable."),
) -> None:
    """
    Run the complete WES analysis pipeline.

    \b
    On first run, ExomeFlow sets itself up automatically:
      • Detects bundled GATK and ANNOVAR
      • Installs missing tools (bwa, samtools, fastp, perl)
      • Asks for reference genome paths (or downloads them)
      • Asks for ANNOVAR database paths (or downloads them)
      • Saves everything — future runs need no extra arguments

    \b
    Workflow
    --------
    FASTQ → fastp → BWA MEM → BAM processing → BQSR
         → HaplotypeCaller → Filtering → ANNOVAR annotation

    \b
    Example
    -------
    exomeflow run --input-dir fastq/ --output results/
    """
    from exomeflow.config import Config
    from exomeflow.pipeline import run_pipeline
    from exomeflow.setup_env import detect_gatk_bin

    # ── First-run setup (skipped if config already complete) ─────────────
    saved = _ensure_ready()

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
        help="ANNOVAR directory. Auto-detected from ExomeFlow folder if omitted.",
    ),
    annovar_db: Optional[Path] = typer.Option(
        None, "--annovar-db",
        help="ANNOVAR humandb directory. Defaults to <annovar-bin>/humandb.",
    ),
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
    from exomeflow.setup_env import detect_annovar_bin, run_setup

    resolved_bin = annovar_bin or detect_annovar_bin()
    if resolved_bin is None:
        console.print(
            "  [red]✘[/red]  ANNOVAR not found. Place the [cyan]annovar/[/cyan] folder "
            "inside the ExomeFlow directory.\n"
            "  Register + download: https://annovar.openbioinformatics.org"
        )
        raise typer.Exit(code=1)

    resolved_db = annovar_db or (resolved_bin / "humandb")

    failed = run_setup(
        refs_dir=refs_dir,
        annovar_bin=resolved_bin,
        annovar_db=resolved_db,
        existing_refs_dir=existing_refs,
    )
    raise typer.Exit(code=min(failed, 1))
