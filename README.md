<div align="center">

<img src="https://raw.githubusercontent.com/imrobintomar/exomeflow/main/ExomeFlow_Icon.png" alt="ExomeFlow Logo" width="700">


# ExomeFlow

**A production-quality Python pipeline for Whole Exome Sequencing analysis.**


| | |
|---|---|
| **CI** | [![CI](https://github.com/imrobintomar/exomeflow/actions/workflows/ci.yml/badge.svg)](https://github.com/imrobintomar/exomeflow/actions/workflows/ci.yml) [![Platform](https://img.shields.io/badge/platform-linux--64-lightgrey)](https://pypi.org/project/exomeflow/) |
| **Package** | [![PyPI Latest Release](https://img.shields.io/pypi/v/exomeflow.svg)](https://pypi.org/project/exomeflow/) [![Downloads](https://static.pepy.tech/badge/exomeflow)](https://pepy.tech/project/exomeflow) [![Wheel](https://img.shields.io/pypi/wheel/exomeflow)](https://pypi.org/project/exomeflow/) [![PyPI Status](https://img.shields.io/pypi/status/exomeflow)](https://pypi.org/project/exomeflow/) |
| **Container** | [![Docker Pulls](https://img.shields.io/docker/pulls/itsrobintomar/exomeflow)](https://hub.docker.com/r/itsrobintomar/exomeflow) [![Docker Image Size](https://img.shields.io/docker/image-size/itsrobintomar/exomeflow/latest)](https://hub.docker.com/r/itsrobintomar/exomeflow) |
| **Meta** | [![License - MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://pypi.org/project/exomeflow/) [![Python Versions](https://img.shields.io/pypi/pyversions/exomeflow)](https://pypi.org/project/exomeflow/) [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.20155767-blue)](https://doi.org/10.5281/zenodo.20155767) [![bio.tools](https://img.shields.io/badge/bio.tools-exomeflow-008080)](https://bio.tools/exomeflow) [![OpenEBench](https://img.shields.io/badge/OpenEBench-benchmarked-blue)](https://openebench.bsc.es/tool/biotools:exomeflow) |
| **Author** | [![AIIMS New Delhi](https://img.shields.io/badge/Institution-AIIMS%20New%20Delhi-red)](https://www.aiims.edu) [![ORCID](https://img.shields.io/badge/ORCID-0009--0002--9084--2019-brightgreen?logo=orcid)](https://orcid.org/0009-0002-9084-2019) |


---

</div>

## Overview

ExomeFlow takes raw paired-end FASTQ files to functionally annotated, clinically
contextualized variants in a single reproducible command. It implements the GATK
best-practices germline short-variant workflow quality control, alignment, base
quality recalibration, variant calling, hard filtering, and ANNOVAR functional
annotation as a modular Python package, with cohort joint genotyping, tumor-only
somatic calling, read-depth CNV detection, and HPO/ACMG clinical annotation available
as opt-in extensions on top of the same pipeline.

The tool is designed to run unattended: on first invocation it detects or provisions
every external dependency it needs (GATK, ANNOVAR, reference genomes, annotation
databases) and persists the resolved configuration, so a single `exomeflow run`
command is sufficient for both a first-time setup and every run after it.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Pipeline Workflow](#pipeline-workflow)
- [Benchmarks](#benchmarks)
- [Installation](#installation)
- [System Requirements](#system-requirements)
- [Python Dependencies](#python-dependencies)
- [Quick Start](#quick-start)
- [Commands](#commands)
- [Reference Files](#reference-files)
- [Input Convention](#input-convention)
- [Output Files](#output-files)
- [Known Limitations](#known-limitations)
- [Getting Help](#getting-help)
- [License](#license)
- [Citation](#citation)

---

## Features

### Core workflow

- **Zero-config first run:**  auto-detects or provisions GATK, ANNOVAR, reference
  genomes, and annotation databases on first use; the resolved configuration is saved
  to `~/.exomeflow/config.json` so every later run needs only `--input-dir`/`--output`
- **GATK best-practice germline calling:** fastp QC → BWA-MEM alignment →
  coordinate sorting → duplicate marking → BQSR → HaplotypeCaller → hard filtering →
  ANNOVAR functional annotation
- **Automatic sample detection:** scans an input directory for paired-end FASTQ
  files by naming convention; no manifest required
- **Per-sample output by default:** any number of samples processed together still
  produces one separate annotated file per sample, identical to running each alone
- **ANNOVAR functional annotation:** against refGene, ClinVar, gnomAD, dbNSFP, ExAC,
  and avSNP150, with automatic per-build (hg38/GRCh37) database selection
- **HPO and ACMG clinical enrichment:** every annotated table is automatically
  joined with HPO gene-to-phenotype terms and ACMG/AMP pathogenicity classification
  (via InterVar); enrichment degrades gracefully rather than blocking a run when its
  dependencies aren't yet available

### Cohort and advanced modes (opt-in)

- **Cohort joint genotyping:** (`--joint-genotyping`) GenomicsDBImport +
  GenotypeGVCFs across all samples, producing one shared cohort VCF and annotation set
- **Somatic tumor-only calling:** (`--mode somatic`) Mutect2 with contamination-aware
  filtering (tumor-normal pairing is not yet supported)
- **Read-depth CNV calling:** (`--cnv`) GATK CollectReadCounts/DenoiseReadCounts/
  PlotDenoisedCopyRatios, no panel-of-normals required
- **GRCh37 or hg38:** (`--genome-build`) reference bundle, ANNOVAR buildver, and
  database selection all follow the requested build automatically
- **Cohort-wide MultiQC rollup:** aggregating fastp, alignment, and GATK metrics,
  generated automatically at the end of every run

### Operational reliability

- **Checkpointed and resumable:** every completed step is recorded per sample; an
  interrupted run resumes exactly where it left off, and upgrading to a version with
  new pipeline steps reprocesses only what's new
- **Parallel cohort processing:** via `--max-workers`
- **Structured logging:** per-sample and pipeline-wide logs with INFO / WARNING /
  ERROR / SUCCESS levels
- **Automatic dependency management:** every external tool and database is
  auto-detected and, if missing, installed or downloaded without a separate setup step
- **Modular architecture:** pipeline steps compose from a pluggable registry, so new
  modes gate in and out without touching unrelated code

---

## Pipeline Workflow

<details>
<summary>Text diagram</summary>

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
                    └──────┬──────┘
                           │  MergeVcfs
                           ▼
┌─────────────────────────────────────────────────────────┐
│  Step 9   SelectVariants  Extract PASS-only variants     │
└──────────────────────────┬──────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────┐
│  Step 10  ANNOVAR         Functional annotation          │
│           refGene · ClinVar · gnomAD · dbNSFP · ExAC     │
│           → multianno.vcf  +  multianno.txt              │
└──────────────────────────┬──────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────┐
│  Step 11  HPO + ACMG      Clinical annotation enrichment │
│           → multianno.hpo.txt                            │
└─────────────────────────────────────────────────────────┘
```

`--joint-genotyping`, `--mode somatic`, and `--cnv` each modify or extend this chain
for the samples they apply to see [Commands](#commands).

</details>

---

## Benchmarks

Benchmarked on **NA12878 (HG001)** whole-exome sequencing data (Hg38), default germline single-sample mode. Accuracy evaluated against the GIAB NISTv4.2.1 .

### Performance

| Metric | Value |
|--------|-------|
| Total runtime | 218.4 min |
| Slowest step | BQSR (141.3 min) |
| Threads | 24 |

### Variant Quality (PASS variants)

| Metric | Value | Expected range |
|--------|-------|----------------|
| SNPs called | 38,413 |  |
| INDELs called | 5,971 |  |
| Ts/Tv ratio | 2.58 | 2.0–3.3 ✓ |
| Het/Hom ratio | 3.10 | 1.5–2.5 |
| dbSNP concordance | 44.7% |  |

### Accuracy (vs. GIAB NISTv4.2.1, PASS-only)

| Variant type | Precision | Recall | F1 score | TP | FP | FN |
|---|---|---|---|---|---|---|
| SNP | 99.41% | 64.67% | 78.36% | 7,787 | 46 | 4,255 |
| INDEL | 89.38% | 66.14% | 76.02% | 623 | 74 | 319 |

> Recall reflects PASS-only evaluation under conservative hard filters. PASS-only
> extraction is unconditional in ExomeFlow  the raw pre-filter VCF is always retained
> alongside it. **This benchmark predates a fix that removed a non-standard
> `DP < 10` site filter** — a threshold that was never part of GATK's own hard-filter
> recommendation and could discard real variants at modest-but-genuine depth (e.g.
> near capture-kit target edges). The recall above likely understates current
> behavior; numbers will be refreshed against a re-run once available.

### Functional Annotation (NA12878)

| Category | Count |
|---|---|
| Total annotated variants | 44,673 |
| Exonic | 15,466 (34.6%) |
| Nonsynonymous SNV | 6,957 |
| Synonymous SNV | 8,158 |
| Stopgain | 57 |
| Frameshift indel | 224 |
| Splicing | 62 |
| ClinVar pathogenic/likely-pathogenic | 5 |
| Novel (not in dbSNP avSNP150) | 658 |

---

## Installation

### Option 1  Python package (recommended)

```bash
pip install exomeflow
```

### Option 2 Docker

```bash
docker pull itsrobintomar/exomeflow:2.2.5

docker run --rm -it \
  -v /path/to/fastq:/data/fastq \
  -v /path/to/refs:/refs \
  -v /path/to/vcf:/vcf \
  -v /path/to/annovar:/annovar \
  -v /path/to/results:/data/results \
  itsrobintomar/exomeflow:2.2.5 run \
    --input-dir    /data/fastq \
    --output       /data/results \
    --reference    /refs/Homo_sapiens_assembly38.fasta \
    --dbsnp        /vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    --mills        /vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-indels /vcf/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --annovar-bin  /annovar \
    --annovar-db   /annovar/humandb \
    --threads      24
```

| Volume mount | Host path | Container path |
|---|---|---|
| Input FASTQs | `/your/fastq/` | `/data/fastq` |
| Reference FASTA + BWA index | `/your/refs/` | `/refs` |
| VCF files (dbSNP, Mills, known indels) | `/your/vcf/` | `/vcf` |
| ANNOVAR scripts | `/your/annovar/` | `/annovar` |
| ANNOVAR humandb | `/your/annovar/humandb/` | `/annovar/humandb` |
| Output | `/your/results/` | `/data/results` |

> ANNOVAR must be mounted  it cannot be bundled due to its licensing terms. Register
> and download at [annovar.openbioinformatics.org](https://annovar.openbioinformatics.org).

### Option 3  Singularity (HPC clusters)

```bash
# Pull directly from Docker Hub
singularity pull exomeflow-2.2.5.sif docker://itsrobintomar/exomeflow:2.2.5

# Or build from the definition file
singularity build exomeflow-2.2.5.sif exomeflow.def

singularity exec \
  --bind /path/to/fastq:/data/fastq \
  --bind /path/to/refs:/refs \
  --bind /path/to/vcf:/vcf \
  --bind /path/to/annovar:/annovar \
  --bind /path/to/results:/data/results \
  exomeflow-2.2.5.sif exomeflow run \
    --input-dir    /data/fastq \
    --output       /data/results \
    --reference    /refs/Homo_sapiens_assembly38.fasta \
    --dbsnp        /vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    --mills        /vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-indels /vcf/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --annovar-bin  /annovar \
    --annovar-db   /annovar/humandb \
    --threads      24
```

<details>
<summary>SLURM job script example</summary>

```bash
#!/bin/bash
#SBATCH --job-name=exomeflow
#SBATCH --cpus-per-task=24
#SBATCH --mem=90G
#SBATCH --time=24:00:00
#SBATCH --output=exomeflow_%j.log

singularity exec \
  --bind $FASTQ_DIR:/data/fastq \
  --bind $REFS_DIR:/refs \
  --bind $VCF_DIR:/vcf \
  --bind $ANNOVAR_DIR:/annovar \
  --bind $RESULTS_DIR:/data/results \
  exomeflow-2.2.5.sif exomeflow run \
    --input-dir    /data/fastq \
    --output       /data/results \
    --reference    /refs/Homo_sapiens_assembly38.fasta \
    --dbsnp        /vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    --mills        /vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-indels /vcf/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --annovar-bin  /annovar \
    --annovar-db   /annovar/humandb \
    --threads      $SLURM_CPUS_PER_TASK
```

</details>

---

## System Requirements

ExomeFlow invokes the following external tools via the command line. Each is
auto-detected on first run and installed automatically where possible (see
[Quick Start](#quick-start)) the table below documents the minimum versions, not a
manual installation requirement.

| Tool | Minimum version | Install |
|------|----------------|---------|
| [BWA](https://github.com/lh3/bwa) | ≥ 0.7.17 | `conda install -c bioconda bwa` |
| [SAMtools](http://www.htslib.org) | ≥ 1.13 | `conda install -c bioconda samtools` |
| [GATK](https://github.com/broadinstitute/gatk/releases) | ≥ 4.6.0 | Auto-downloaded, or `conda install -c bioconda gatk4` |
| [fastp](https://github.com/OpenGENOMICS/fastp) | ≥ 0.20.1 | `conda install -c bioconda fastp` |
| [Perl](https://www.perl.org) | ≥ 5.26 | `conda install perl` |
| [ANNOVAR](https://annovar.openbioinformatics.org) | latest | Requires registration; not auto-installable |

> ANNOVAR is the one dependency ExomeFlow cannot provision for you  its license
> requires individual registration. Everything else, including GATK itself, is
> detected or fetched automatically.

---

## Python Dependencies

- **[typer](https://typer.tiangolo.com/)** CLI interface
- **[rich](https://rich.readthedocs.io/)** Structured, coloured terminal output
- **[pandas](https://pandas.pydata.org/)** Variant summaries and HPO/ACMG annotation joins

Installed automatically with `pip install exomeflow`. `matplotlib` (needed only for
`--cnv` plots) is an optional extra (`pip install exomeflow[viz]`) that the dependency
checker installs automatically the first time `--cnv` is used.

---

## Quick Start

```bash
pip install exomeflow
```

Arrange paired-end FASTQ files following the [naming convention](#input-convention):

```
fastq/
├── sample1_1.fastq.gz
├── sample1_2.fastq.gz
├── sample2_1.fastq.gz
└── sample2_2.fastq.gz
```

Run the pipeline:

```bash
exomeflow run --input-dir fastq/ --output results/
```

On first run, ExomeFlow detects or provisions GATK, ANNOVAR, reference data, the HPO
gene-to-phenotype mapping, and InterVar, then saves the resolved configuration to
`~/.exomeflow/config.json`. Every later run needs nothing beyond `--input-dir` and
`--output`. Add `--yes` to skip all interactive prompts for unattended or CI use.

To control every path explicitly instead for example on a shared cluster where
reference data already exists at known locations:

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

`exomeflow setup` is also available to run provisioning as its own step, independent
of a pipeline run useful for switching reference builds or refreshing databases.

### Cohort, somatic, CNV, and GRCh37 modes

```bash
# Cohort joint genotyping instead of per-sample VCFs
exomeflow run --input-dir fastq/ --output results/ --joint-genotyping --intervals targets.bed

# Somatic tumor-only calling with Mutect2
exomeflow run --input-dir fastq/ --output results/ --mode somatic \
  --germline-resource af-only-gnomad.vcf.gz --panel-of-normals pon.vcf.gz

# Read-depth CNV calling alongside the standard germline workflow
exomeflow run --input-dir fastq/ --output results/ --cnv --intervals targets.bed

# GRCh37/hg19 instead of hg38
exomeflow run --input-dir fastq/ --output results/ --genome-build GRCh37
```

---

## Commands

### `exomeflow run` Execute the WES pipeline

```
exomeflow run [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--input-dir`, `-i` | required | Directory containing paired FASTQ files |
| `--output`, `-o` | `results/` | Root output directory |
| `--reference`, `-r` | auto-resolved | BWA-indexed reference FASTA |
| `--dbsnp` | auto-resolved | dbSNP VCF (bgzipped + tabix-indexed) |
| `--mills` | auto-resolved | Mills and 1000G gold standard indels VCF |
| `--known-indels` | auto-resolved | Known indels VCF for BQSR |
| `--intervals` | optional | Exome capture BED file — required for `--joint-genotyping`/`--cnv` |
| `--interval-padding` | `100` | Base-pair padding around each target interval |
| `--annovar-bin` | auto-resolved | Directory containing `table_annovar.pl` |
| `--annovar-db` | auto-resolved | ANNOVAR humandb directory |
| `--annovar-protocols` | build default | Override the ANNOVAR `--protocol` list |
| `--annovar-operations` | build default | Override the ANNOVAR `--operation` list (must pair with `--annovar-protocols`) |
| `--mode` | `germline` | `germline` (HaplotypeCaller) or `somatic` (tumor-only Mutect2) |
| `--genome-build` | `hg38` | `hg38` or `GRCh37` |
| `--joint-genotyping` | off | Cohort mode: one shared VCF/annotation instead of per-sample files |
| `--cnv` | off | Also call read-depth CNVs per sample (needs `--intervals`) |
| `--germline-resource` | optional | gnomAD AF-only VCF for Mutect2, used by `--mode somatic` |
| `--panel-of-normals` | optional | Pre-built Panel of Normals VCF for Mutect2, used by `--mode somatic` |
| `--threads`, `-t` | auto-detected | Threads for BWA MEM and GATK HaplotypeCaller |
| `--fastp-threads` | `8` | Threads for fastp |
| `--annovar-threads` | `24` | Threads for ANNOVAR |
| `--max-workers` | `1` | Number of samples to process in parallel |
| `--java-opts` | auto-sized | JVM options passed via `JAVA_OPTS` |
| `--yes`, `-y` | off | Skip all interactive prompts (unattended/CI runs) |

`--threads` and `--java-opts` are sized automatically from available CPU/RAM when not
set explicitly.

### `exomeflow setup` Run provisioning independently of a pipeline run

```
exomeflow setup [--refs-dir PATH] [--annovar-bin PATH] [--annovar-db PATH] [--genome-build hg38|GRCh37] [--existing-refs PATH] [--yes]
```

Not required before `exomeflow run` first-run auto-setup covers the same ground.
Useful for re-provisioning without running the pipeline itself.

---

## Reference Files

| File | Source | Size |
|------|--------|------|
| Reference FASTA + BWA index | GATK resource bundle | ~10 GB |
| dbSNP VCF | GATK resource bundle | ~10 GB |
| Mills and 1000G gold standard indels | GATK resource bundle | ~200 MB |
| Known indels VCF | GATK resource bundle | ~100 MB |
| Exome capture BED | Your sequencing kit vendor | varies |
| ANNOVAR humandb | ANNOVAR download server | ~90 GB |

`exomeflow run`/`exomeflow setup` download all GATK resource bundle files
automatically for the requested `--genome-build`. Manual download, if needed:

```bash
gsutil -m cp -r gs://gcp-public-data--broad-references/hg38/v0/ /data/refs/
```

---

## Input Convention

ExomeFlow detects samples automatically from paired-end FASTQ filenames:

```
<sample_id>_1.fastq.gz   ← Read 1
<sample_id>_2.fastq.gz   ← Read 2
```

`sample_id` can be any string sharing a common prefix between the two read files 
an SRR accession, a patient ID, or anything else.

---

## Output Files

Per-sample output (default one complete set per sample, regardless of how many
samples are in the run):

| File | Description |
|------|-------------|
| `Mapsam/<sample>_recalibrated.bam` | Analysis-ready BAM open in IGV |
| `VCF/<sample>.vcf` | Raw HaplotypeCaller output (germline) |
| `VCF/<sample>_unfiltered.vcf.gz` | Raw Mutect2 output (`--mode somatic`) |
| `VCF/<sample>_PASS.vcf` | PASS-only filtered variants |
| `VCF/<sample>.annovar.<buildver>_multianno.{vcf,txt}` | ANNOVAR-annotated variants |
| `VCF/<sample>.annovar.hpo.txt` | Annotated table with HPO terms and ACMG classification |
| `filtered_fastp/<sample>_fastp.html` | fastp QC report |
| `Mapsam/<sample>_flagstat.txt` | Alignment statistics |
| `CNV/<sample>_denoised_cr.tsv` + plot | Read-depth CNV calls (`--cnv` only) |
| `logs/analysis_<timestamp>.log` | Full pipeline log |
| `logs/<sample>_<timestamp>.log` | Per-sample log |

Cohort output (`--joint-genotyping` only replaces the per-sample VCF/annotation
files above with one shared set):

| File | Description |
|------|-------------|
| `VCF/cohort/cohort.vcf.gz` | Joint-genotyped multi-sample VCF |
| `VCF/cohort/cohort_PASS.vcf` | PASS-only filtered cohort variants |
| `VCF/cohort/cohort.annovar.<buildver>_multianno.{vcf,txt}` | Annotated cohort variants |
| `VCF/cohort/cohort.annovar.hpo.txt` | Annotated cohort table with HPO/ACMG |

Generated at the end of every run:

| File | Description |
|------|-------------|
| `multiqc/exomeflow_report.html` | Cohort-wide QC rollup (fastp, flagstat, GATK metrics) |

---

## Known Limitations

- **Somatic mode is tumor-only.** Matched tumor-normal pairing is not yet supported.
- **HPO/ACMG enrichment depends on external, unversioned resources** (the HPO
  gene-to-phenotype mapping and InterVar's own reference databases). When these
  aren't available, the corresponding columns are omitted rather than the run failing
  — check the pipeline log if `HPO_terms`/`ACMG_classification` are missing from a
  given sample's output.
- **GRCh37/hg19 support annotates with gnomAD v2.1.1**, not v4.1 gnomAD v4 was never
  released for that build.

---

---

## License

MIT — see [pypi.org/project/exomeflow](https://pypi.org/project/exomeflow/) for details.

---
## Getting Help

For usage questions and bug reports, contact:

**Robin Kumar**  itsrobintomar@gmail.com



---
**Made with ❤️ ||  Dr Prabudh Goel lab, AIIMS New Delhi , India**
