# ExomeFlow

**A production-quality Python pipeline for Whole Exome Sequencing analysis.**

ExomeFlow takes raw paired-end FASTQ files to functionally annotated, clinically
contextualized variants in a single reproducible command. It implements the GATK
best-practices germline short-variant workflow — quality control, alignment, base
quality recalibration, variant calling, hard filtering, and ANNOVAR functional
annotation — as a modular Python package, with cohort joint genotyping, tumor-only
somatic calling, read-depth CNV detection, and HPO/ACMG clinical annotation available
as opt-in extensions on top of the same pipeline.

The tool is designed to run unattended: on first invocation it detects or provisions
every external dependency it needs (GATK, ANNOVAR, reference genomes, annotation
databases) and persists the resolved configuration, so a single `exomeflow run`
command is sufficient for both a first-time setup and every run after it.

```bash
pip install exomeflow
exomeflow run --input-dir fastq/ --output results/
```

## Where to go next

- **[Installation](installation.md)** — install the package and let it provision
  everything else on first run.
- **[CLI Reference](usage.md)** — every flag `exomeflow run` accepts.
- **[Modes & Options](modes.md)** — cohort joint genotyping, somatic calling,
  CNV detection, GRCh37 support.
- **[Architecture](architecture.md)** — how the pluggable step registry composes
  a pipeline from `Config` flags.
- **[Real-World Example — NA12878](example-na12878.md)** — a full run against a
  public reference sample, real commands, real output.
- **[Output Files](output-files.md)** — what each pipeline step writes and where.
- **[Troubleshooting](troubleshooting.md)** — common errors and fixes.

## Project links

- [GitHub repository](https://github.com/imrobintomar/exomeflow)
- [PyPI package](https://pypi.org/project/exomeflow/)
- [Changelog](changelog.md)

## Known limitations

- Somatic mode is tumor-only; tumor-normal paired calling is not yet supported.
- HPO term annotation and ACMG classification (InterVar) depend on external
  databases/tools that are auto-provisioned on first run — if that provisioning
  fails (e.g. no network), those two columns are skipped with a warning rather
  than failing the run.
- GRCh37 annotation uses gnomAD v2.1.1 (gnomAD v4.1 was never released for hg19).
