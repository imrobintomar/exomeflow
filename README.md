![ExomeFlow Cover](https://raw.githubusercontent.com/imrobintomar/exomeflow/main/ExomeFlow_Cover.png)

# ExomeFlow: A Production-Quality Python WES Analysis Toolkit

![ExomeFlow Icon](https://raw.githubusercontent.com/imrobintomar/exomeflow/main/ExomeFlow_Icon.png)

| | |
|---|---|
| **Testing** | [![CI](https://img.shields.io/badge/CI-passing-brightgreen)](https://github.com/imrobintomar/exomeflow/actions) |
| **Package** | [![PyPI Latest Release](https://img.shields.io/pypi/v/exomeflow.svg)](https://pypi.org/project/exomeflow/) [![PyPI Downloads](https://img.shields.io/pypi/dm/exomeflow)](https://pypi.org/project/exomeflow/) |
| **Meta** | [![License - MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/imrobintomar/exomeflow/blob/main/LICENSE) [![Python Versions](https://img.shields.io/pypi/pyversions/exomeflow)](https://pypi.org/project/exomeflow/) [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo-blue)](https://github.com/imrobintomar/exomeflow) |

---

## What is it?

**ExomeFlow** is a Python package that provides a complete, automated Whole Exome Sequencing (WES)
analysis workflow from raw FASTQ files to functionally annotated variants in a single
reproducible CLI command.

It aims to be the standard high-level pipeline for WES analysis in Python, combining
GATK best-practice variant calling, hard filtering, and ANNOVAR annotation into one
modular, maintainable package. It handles cohort-level processing (multiple samples),
checkpointing for resumable runs, structured logging, and parallel execution out of the box.

---

## Table of Contents

- [What is it?](#what-is-it)
- [Main Features](#main-features)
- [Pipeline Workflow](#pipeline-workflow)
- [Where to get it](#where-to-get-it)
- [System Requirements](#system-requirements)
- [Python Dependencies](#python-dependencies)
- [Quick Start](#quick-start)
- [Commands](#commands)
- [Reference Files](#reference-files)
- [Input Convention](#input-convention)
- [Output Files](#output-files)
- [Documentation](#documentation)
- [Getting Help](#getting-help)
- [License](#license)
- [Citation](#citation)

---

## Main Features

Here are the things ExomeFlow does well:

- **One-command setup** — `exomeflow setup` installs all system tools, downloads hg38 reference
  files (~13 GB) and ANNOVAR databases (~100 GB) automatically
- **Automatic sample detection** — scans an input directory and detects all paired-end
  samples from FASTQ filenames; no manifest file required
- **Complete GATK best-practice workflow** — fastp QC → BWA MEM alignment → coordinate
  sorting → duplicate marking → BQSR → HaplotypeCaller → hard filtering → ANNOVAR annotation
- **Cohort processing** — processes any number of samples sequentially or in parallel
  with `--max-workers`
- **Checkpointing and resume** — every completed step is recorded; an interrupted run
  resumes exactly where it left off without repeating work
- **Automatic requirements check** — verifies all system tools and Python packages
  before the pipeline starts, reporting every missing dependency at once
- **Structured logging** — per-sample log files plus a pipeline-wide log with
  INFO / WARNING / ERROR / SUCCESS levels
- **GATK hard filters** — applies GATK best-practice SNP and INDEL hard-filter
  thresholds and extracts PASS-only variants automatically
- **ANNOVAR functional annotation** — annotates variants against 8 databases:
  refGene, ClinVar, gnomAD, dbNSFP, COSMIC, ExAC, avSNP150, and dbscSNV
- **Modular architecture** — each pipeline step is an independent Python module;
  easy to extend or modify individual steps without touching the rest
- **PyPI installable** — `pip install exomeflow`; no Docker or Nextflow required

---

## Pipeline Workflow

```
Raw FASTQ
    │
    ▼
┌─────────────────────────────────────────────────────────┐
│  Step 1   fastp         Quality control & adapter trim   │
│           length ≥ 50 bp · base quality ≥ Q30            │
└──────────────────────────┬──────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────┐
│  Step 2   BWA MEM        Read alignment to hg38          │
│           -Y -K 100000000 · read-group tags set          │
└──────────────────────────┬──────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────┐
│  Step 3   GATK SortSam   Coordinate-sort BAM             │
│  Step 4   samtools       Flagstat alignment QC           │
│  Step 5   GATK MarkDuplicates   PCR duplicate removal    │
│  Step 6   GATK BuildBamIndex    BAI index                │
└──────────────────────────┬──────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────┐
│  Step 7   GATK BQSR      BaseRecalibrator + ApplyBQSR    │
│           Known sites: dbSNP · Mills · known indels      │
│           → recalibrated.bam  (IGV-ready)                │
└──────────────────────────┬──────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────┐
│  Step 8   GATK HaplotypeCaller   Variant calling         │
│           Exome intervals + padding · dbSNP annotation   │
└──────────────────────────┬──────────────────────────────┘
                           │
                    ┌──────┴──────┐
                    ▼             ▼
               SNP filters   INDEL filters
               (Step 9)       (Step 10)
                    └──────┬──────┘
                           │  MergeVcfs
                           ▼
┌─────────────────────────────────────────────────────────┐
│  Step 11  SelectVariants  Extract PASS-only variants     │
└──────────────────────────┬──────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────┐
│  Step 12  ANNOVAR         Functional annotation          │
│           refGene · ClinVar · gnomAD · dbNSFP · COSMIC   │
│           → multianno.vcf  +  multianno.txt              │
└─────────────────────────────────────────────────────────┘
```

---

## Where to get it

The source code is hosted on GitHub at:
**https://github.com/imrobintomar/exomeflow**

ExomeFlow is available via three installation methods:

### Option 1 — Python Package (recommended)

```bash
pip install exomeflow
```

### Option 2 — Docker

```bash
# Pull image
docker pull imrobintomar/exomeflow:latest

# Run pipeline
docker run --rm -it \
  -v /path/to/fastq:/data/fastq \
  -v /path/to/refs:/refs \
  -v /path/to/annovar:/annovar \
  -v /path/to/results:/data/results \
  imrobintomar/exomeflow:latest run \
    --input-dir    /data/fastq \
    --output       /data/results \
    --reference    /refs/hg38.fa \
    --dbsnp        /refs/dbsnp.vcf.gz \
    --mills        /refs/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-indels /refs/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --annovar-bin  /annovar \
    --annovar-db   /annovar/humandb \
    --threads      24
```

### Option 3 — Singularity (HPC clusters)

```bash
# Pull from Docker Hub
singularity pull docker://imrobintomar/exomeflow:latest

# Run pipeline
singularity exec exomeflow_latest.sif exomeflow run \
  --input-dir    /path/to/fastq \
  --output       /path/to/results \
  --reference    /path/to/hg38.fa \
  --dbsnp        /path/to/dbsnp.vcf.gz \
  --mills        /path/to/mills.vcf.gz \
  --known-indels /path/to/known_indels.vcf.gz \
  --annovar-bin  /path/to/annovar \
  --annovar-db   /path/to/annovar/humandb \
  --threads      24
```

> **Note:** ANNOVAR requires registration at [annovar.openbioinformatics.org](https://annovar.openbioinformatics.org)
> and must be mounted as a volume (`-v /your/annovar:/annovar`).
> It cannot be bundled in the Docker image due to licensing restrictions.

The list of changes between each release can be found in the
[Release History](https://github.com/imrobintomar/exomeflow/releases).

---

## System Requirements

ExomeFlow calls the following external tools via the command line.
They must be installed separately and available on your `PATH`.

| Tool | Minimum Version | Install |
|------|----------------|---------|
| [BWA](https://github.com/lh3/bwa) | ≥ 0.7.17 | `conda install -c bioconda bwa` |
| [SAMtools](http://www.htslib.org) | ≥ 1.13 | `conda install -c bioconda samtools` |
| [GATK](https://github.com/broadinstitute/gatk/releases) | ≥ 4.6.0 | `conda install -c bioconda gatk4` |
| [fastp](https://github.com/OpenGENOMICS/fastp) | ≥ 0.20.1 | `conda install -c bioconda fastp` |
| [Perl](https://www.perl.org) | ≥ 5.26 | `conda install perl` |
| [ANNOVAR](https://annovar.openbioinformatics.org) | latest | Register + download from website |

> **Tip:** Run `exomeflow setup` after installation to automatically verify tools,
> download hg38 reference files, and populate ANNOVAR databases in one step.

---

## Python Dependencies

- **[typer](https://typer.tiangolo.com/)** — Builds the CLI interface
- **[rich](https://rich.readthedocs.io/)** — Provides coloured terminal output and structured logging
- **[pandas](https://pandas.pydata.org/)** — Data handling for variant count summaries

All Python dependencies are installed automatically with `pip install exomeflow`.

---

## Quick Start

### 1. Install ExomeFlow

```bash
pip install exomeflow
```

### 2. Set up all dependencies and reference data

```bash
exomeflow setup \
  --refs-dir   /data/references/hg38 \
  --annovar-bin /opt/annovar \
  --annovar-db  /opt/annovar/humandb
```

This command will:
- Install missing Python packages
- Install system tools via conda (fastp, bwa, samtools, gatk4, perl)
- Download hg38 reference files (~13 GB) using gsutil or wget
- Download ANNOVAR annotation databases (~100 GB)

### 3. Prepare FASTQ files

```
fastq/
├── sample1_1.fastq.gz
├── sample1_2.fastq.gz
├── sample2_1.fastq.gz
└── sample2_2.fastq.gz
```

### 4. Run the pipeline

```bash
exomeflow run \
  --input-dir    fastq/ \
  --output       results/ \
  --reference    /data/references/hg38/hg38.fa \
  --dbsnp        /data/references/hg38/dbsnp.vcf.gz \
  --mills        /data/references/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known-indels /data/references/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --intervals    refs/Illumina_Exome_TargetedRegions_v1.2.hg38.bed \
  --annovar-bin  /opt/annovar \
  --annovar-db   /opt/annovar/humandb \
  --threads      32 \
  --max-workers  2
```

---

## Commands

### `exomeflow setup` — Install dependencies and download reference data

```
exomeflow setup --refs-dir PATH --annovar-bin PATH --annovar-db PATH
```

| Option | Description |
|--------|-------------|
| `--refs-dir` | Directory to download hg38 reference files into |
| `--annovar-bin` | ANNOVAR installation directory (must contain `annotate_variation.pl`) |
| `--annovar-db` | ANNOVAR humandb directory for database downloads |

### `exomeflow run` — Execute the WES pipeline

```
exomeflow run [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--input-dir`, `-i` | required | Directory containing paired FASTQ files |
| `--output`, `-o` | `results/` | Root output directory |
| `--reference`, `-r` | required | BWA-indexed reference FASTA (hg38.fa) |
| `--dbsnp` | required | dbSNP VCF (bgzipped + tabix-indexed) |
| `--mills` | required | Mills and 1000G gold standard indels VCF |
| `--known-indels` | required | Known indels VCF for BQSR |
| `--intervals` | _(optional)_ | Exome capture BED file |
| `--interval-padding` | `100` | Base-pair padding around each target interval |
| `--annovar-bin` | required | Directory containing `table_annovar.pl` |
| `--annovar-db` | required | ANNOVAR humandb directory |
| `--threads`, `-t` | `24` | Threads for BWA MEM and GATK HaplotypeCaller |
| `--fastp-threads` | `8` | Threads for fastp |
| `--annovar-threads` | `24` | Threads for ANNOVAR |
| `--max-workers` | `1` | Number of samples to process in parallel |
| `--java-opts` | `-Xmx80g` | JVM options passed via JAVA_OPTS |

---

## Reference Files

| File | Source | Size |
|------|--------|------|
| `hg38.fa` + BWA index | UCSC / GATK resource bundle | ~10 GB |
| `dbsnp.vcf.gz` | GATK resource bundle | ~10 GB |
| `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` | GATK resource bundle | ~200 MB |
| `Homo_sapiens_assembly38.known_indels.vcf.gz` | GATK resource bundle | ~100 MB |
| Exome capture BED | Your sequencing kit vendor | varies |
| ANNOVAR humandb (8 databases) | ANNOVAR download server | ~100 GB |

`exomeflow setup` downloads all GATK resource bundle files automatically.

Manual download:

```bash
gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38/v0/ /data/refs/
```

---

## Input Convention

ExomeFlow automatically detects samples from paired-end FASTQ filenames.
Files must follow the pattern:

```
<sample_id>_1.fastq.gz   ← Read 1
<sample_id>_2.fastq.gz   ← Read 2
```

The `sample_id` can be any string — SRR accession, patient ID, etc.

---

## Output Files

| File | Description |
|------|-------------|
| `Mapsam/<sample>_recalibrated.bam` | Analysis-ready BAM — open in IGV |
| `VCF/<sample>.vcf` | Raw HaplotypeCaller output |
| `VCF/<sample>_PASS.vcf` | PASS-only hard-filtered variants |
| `VCF/<sample>.annovar.hg38_multianno.vcf` | Annotated VCF |
| `VCF/<sample>.annovar.hg38_multianno.txt` | Annotated tab-delimited table |
| `filtered_fastp/<sample>_fastp.html` | fastp QC report |
| `Mapsam/<sample>_flagstat.txt` | Alignment statistics |
| `logs/analysis_<timestamp>.log` | Full pipeline log |
| `logs/<sample>_<timestamp>.log` | Per-sample log |

---

## Documentation

Full usage documentation is available in [USAGE.md](USAGE.md), including:

- Complete CLI option reference
- How to resume interrupted runs
- How to tune parallel processing
- Common errors and fixes
- Quick reference card

---

## Getting Help

For usage questions, please open a
[GitHub Issue](https://github.com/imrobintomar/exomeflow/issues).

Bug reports, feature requests, and general questions are all welcome.

---

## License

[MIT](LICENSE)

---

## Citation

If you use ExomeFlow in your research, please cite:

> Robin Tomar. *ExomeFlow: A Production-Quality Python Package for Automated
> Whole Exome Sequencing Analysis*. AIIMS New Delhi, 2025.
> https://pypi.org/project/exomeflow/

---

*Built for the bioinformatics community · Robin Tomar, AIIMS New Delhi*
