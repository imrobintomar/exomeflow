"""
Typer-based CLI for ExomeFlow.

Entry point registered in pyproject.toml:
    exomeflow = "exomeflow.cli:app"
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from exomeflow import __version__

app = typer.Typer(
    name="exomeflow",
    help="Whole Exome Sequencing analysis pipeline.",
    add_completion=False,
    rich_markup_mode="rich",
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
    """ExomeFlow — production WES pipeline."""


@app.command("run")
def run_command(
    # ------------------------------------------------------------------ I/O
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
    # ---------------------------------------------------------- reference files
    reference: Path = typer.Option(
        ...,
        "--reference", "-r",
        help="Path to the BWA-indexed reference genome FASTA (e.g. hg38.fa).",
    ),
    dbsnp: Path = typer.Option(
        ...,
        "--dbsnp",
        help="Path to the dbSNP VCF (bgzipped + tabix-indexed).",
    ),
    mills: Path = typer.Option(
        ...,
        "--mills",
        help="Path to Mills and 1000G gold standard indels VCF.",
    ),
    known_indels: Path = typer.Option(
        ...,
        "--known-indels",
        help="Path to Homo sapiens assembly known indels VCF.",
    ),
    intervals: Path = typer.Option(
        Path(""),
        "--intervals",
        help="Exome capture BED file. Omit to run in whole-genome mode.",
    ),
    interval_padding: int = typer.Option(
        100,
        "--interval-padding",
        help="Base-pair padding around each target interval.",
    ),
    # ------------------------------------------------------------ ANNOVAR
    annovar_bin: Path = typer.Option(
        ...,
        "--annovar-bin",
        help="Directory containing table_annovar.pl (ANNOVAR installation).",
    ),
    annovar_db: Path = typer.Option(
        ...,
        "--annovar-db",
        help="ANNOVAR humandb directory for hg38.",
    ),
    # ----------------------------------------------------------- performance
    threads: int = typer.Option(
        24,
        "--threads", "-t",
        help="Threads for BWA MEM and GATK HaplotypeCaller.",
        min=1,
    ),
    fastp_threads: int = typer.Option(
        8,
        "--fastp-threads",
        help="Threads for fastp.",
        min=1,
    ),
    annovar_threads: int = typer.Option(
        24,
        "--annovar-threads",
        help="Threads for ANNOVAR.",
        min=1,
    ),
    max_workers: int = typer.Option(
        1,
        "--max-workers",
        help="Number of samples to process in parallel.",
        min=1,
    ),
    java_opts: str = typer.Option(
        "-Xmx80g",
        "--java-opts",
        help="JVM options passed via JAVA_OPTS environment variable.",
    ),
) -> None:
    """
    Run the complete WES analysis pipeline.

    \b
    Workflow
    --------
    FASTQ → fastp → BWA MEM → BAM processing → BQSR
         → HaplotypeCaller → Filtering → ANNOVAR annotation

    \b
    Example
    -------
    exomeflow run \\
      --input-dir fastq/ \\
      --reference hg38.fa \\
      --dbsnp dbsnp.vcf.gz \\
      --mills mills.vcf.gz \\
      --known-indels known_indels.vcf.gz \\
      --intervals exome_targets.bed \\
      --annovar-db /path/to/annovar/humandb \\
      --annovar-bin /path/to/annovar \\
      --threads 32 \\
      --max-workers 2 \\
      --output results/
    """
    # Late import keeps startup fast
    from exomeflow.config import Config
    from exomeflow.pipeline import run_pipeline

    cfg = Config(
        input_dir=input_dir,
        output_dir=output,
        reference=reference,
        dbsnp=dbsnp,
        mills=mills,
        known_indels=known_indels,
        intervals=intervals,
        interval_padding=interval_padding,
        annovar_bin=annovar_bin,
        annovar_db=annovar_db,
        threads=threads,
        fastp_threads=fastp_threads,
        annovar_threads=annovar_threads,
        max_workers=max_workers,
        java_opts=java_opts,
    )

    failed = run_pipeline(cfg)
    raise typer.Exit(code=min(failed, 1))
