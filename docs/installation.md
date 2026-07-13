# Installation

## 1. Install the package

```bash
pip install exomeflow
```

This installs the `exomeflow` CLI and its Python dependencies (`typer`, `rich`,
`pandas`; `matplotlib` is pulled in automatically only if you use `--cnv`).

## 2. Run it

```bash
exomeflow run --input-dir fastq/ --output results/
```

That's the whole installation story from the user's side. On first run, ExomeFlow
walks through setup automatically:

1. Detects GATK and ANNOVAR on `PATH`, or auto-downloads/clones them if missing.
2. Installs missing system tools (`bwa`, `samtools`, `fastp`, `perl`) via `conda`.
3. Asks for hg38/GRCh37 reference genome paths — or downloads the GATK resource
   bundle for you.
4. Asks for ANNOVAR database paths — or downloads the missing databases.
5. Auto-provisions HPO gene-to-phenotype data, InterVar (ACMG classification),
   and MultiQC.
6. Sizes `--threads` / `--java-opts` from your machine's actual CPU count and RAM.
7. Saves every resolved path to `~/.exomeflow/config.json` — every run after this
   one needs no extra arguments unless you want to override something.

Add `--yes` to auto-confirm every setup prompt instead of asking interactively —
useful for CI or unattended/background runs:

```bash
exomeflow run --input-dir fastq/ --output results/ --yes
```

## Re-running setup

If you need to change a reference path, point at a different ANNOVAR database
directory, or repair a partially-provisioned config, use the `setup` subcommand
directly instead of `run`:

```bash
exomeflow setup
```

## FASTQ naming convention

ExomeFlow discovers samples automatically from paired-end FASTQ files matching
`<sample_id>_1.fastq.gz` / `<sample_id>_2.fastq.gz` in `--input-dir`:

```
fastq/
├── sample1_1.fastq.gz
├── sample1_2.fastq.gz
├── sample2_1.fastq.gz
└── sample2_2.fastq.gz
```

`sample_id` can be any string (SRR accession, patient ID, etc.) as long as both
read files share the same prefix.

## Requirements

- Linux, Python ≥ 3.9.
- Enough disk space for the hg38/GRCh37 reference bundle and ANNOVAR databases
  (tens of GB) if you don't already have them — ExomeFlow downloads what's
  missing on first run.
