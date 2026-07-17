# ExomeFlow — How to Use

**Author:** Robin Kumar, Bioinformatics Scientist, AIIMS New Delhi  
**Package:** `pip install exomeflow`  
**Version:** 2.1.2

---

## Table of Contents

1. [Installation](#1-installation)
2. [Before You Start — Requirements Check](#2-before-you-start--requirements-check)
3. [Prepare Your Input Files](#3-prepare-your-input-files)
4. [Run the Pipeline](#4-run-the-pipeline)
5. [All Command-Line Options](#5-all-command-line-options)
6. [Understanding the Output](#6-understanding-the-output)
7. [Resuming an Interrupted Run](#7-resuming-an-interrupted-run)
8. [Running Multiple Samples in Parallel](#8-running-multiple-samples-in-parallel)
9. [Common Errors and Fixes](#9-common-errors-and-fixes)

---

## 1. Installation

### Install ExomeFlow

```bash
pip install exomeflow
```

### Install system tools (if not already installed)

ExomeFlow calls these tools via the command line.  
They must be installed separately and available on your `PATH`.

```bash
# BWA
conda install -c bioconda bwa

# SAMtools
conda install -c bioconda samtools

# fastp
conda install -c bioconda fastp

# GATK — download from GitHub and add to PATH
# https://github.com/broadinstitute/gatk/releases
export PATH="/path/to/gatk-4.6.2.0:$PATH"

# Perl (usually pre-installed)
conda install perl

# ANNOVAR — register and download from:
# https://annovar.openbioinformatics.org
```

### Verify everything is installed correctly

```bash
python check_requirements.py
```

Expected output:

```
══════════════════════════════════════════════════════════
  ExomeFlow — Requirements Check
══════════════════════════════════════════════════════════

Python Interpreter
  ✔  Python 3.13.x  (required >= 3.9)

Python Packages
  ✔  typer 0.20.0   (required >= 0.12.0)
  ✔  rich  14.2.0   (required >= 13.0.0)
  ✔  pandas 3.0.0   (required >= 2.0.0)
  ✔  matplotlib 3.10.x  (required >= 3.7.0)

System Tools
  ✔  bwa      0.7.17  (required >= 0.7.17)
  ✔  samtools 1.13    (required >= 1.13)
  ✔  gatk     4.6.2   (required >= 4.6.0)
  ✔  fastp    0.20.1  (required >= 0.20.1)
  ✔  perl     5.34.0  (required >= 5.26)
  ✔  ANNOVAR found at /path/to/annovar/table_annovar.pl

══════════════════════════════════════════════════════════
  ✔  All requirements satisfied — ExomeFlow is ready to run.
══════════════════════════════════════════════════════════
```

> **Note:** The requirements check also runs **automatically** as Step 0
> every time you execute `exomeflow run`. If anything is missing, the
> pipeline will stop immediately and tell you exactly what to fix.

---

## 2. Before You Start — Requirements Check

You will need the following reference files for hg38:

| File | Where to get it |
|------|----------------|
| `hg38.fa` | UCSC / Ensembl / GATK resource bundle |
| `dbsnp.vcf.gz` + `.tbi` | GATK hg38 resource bundle |
| `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` | GATK hg38 resource bundle |
| `Homo_sapiens_assembly38.known_indels.vcf.gz` | GATK hg38 resource bundle |
| Exome capture BED | From your capture kit vendor (Illumina / Twist) |
| ANNOVAR humandb | `annotate_variation.pl -buildver hg38 -downdb` |

GATK resource bundle download:

```bash
# Using gsutil (Google Cloud)
gsutil -m cp -r gs://gatk-best-practices/somatic-hg38/ .
```

---

## 3. Prepare Your Input Files

### FASTQ naming convention

ExomeFlow automatically detects samples from paired-end FASTQ files.  
Two naming conventions are accepted — pick whichever matches your data:

```
fastq/
├── sample1_1.fastq.gz      ← Read 1
├── sample1_2.fastq.gz      ← Read 2
├── sample2_1.fastq.gz
├── sample2_2.fastq.gz
├── sample3_R1.fastq.gz     ← Read 1 (alternate convention)
└── sample3_R2.fastq.gz     ← Read 2 (alternate convention)
```

Pattern: `<sample_id>_1.fastq.gz` / `<sample_id>_2.fastq.gz`,
or `<sample_id>_R1.fastq.gz` / `<sample_id>_R2.fastq.gz`.

The `sample_id` can be any string (SRR accession, patient ID, etc.)  
as long as both read files share the same prefix and convention. A single
run can mix samples using different conventions (e.g. `sample1_1.fastq.gz`
alongside `sample2_R1.fastq.gz`) — the convention is resolved per sample.

### Recommended directory layout

```
project/
├── fastq/                          ← Put your FASTQ files here
│   ├── sample1_1.fastq.gz
│   └── sample1_2.fastq.gz
├── refs/                           ← Reference files
│   ├── hg38.fa
│   ├── hg38.fa.fai
│   ├── hg38.dict
│   ├── dbsnp.vcf.gz
│   ├── dbsnp.vcf.gz.tbi
│   ├── Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
│   └── Homo_sapiens_assembly38.known_indels.vcf.gz
├── intervals/
│   └── Illumina_Exome_TargetedRegions_v1.2.hg38.bed
└── results/                        ← ExomeFlow writes output here
```

---

## 4. Run the Pipeline

### Minimal command

On first run, ExomeFlow auto-detects/installs everything it needs and saves the
resolved paths to `~/.exomeflow/config.json` — so this is genuinely all you need:

```bash
exomeflow run --input-dir fastq/ --output results/
```

### Explicit command (override auto-resolved paths)

```bash
exomeflow run \
  --input-dir    fastq/ \
  --output       results/ \
  --reference    refs/hg38.fa \
  --dbsnp        refs/dbsnp.vcf.gz \
  --mills        refs/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known-indels refs/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --intervals    intervals/Illumina_Exome_TargetedRegions_v1.2.hg38.bed \
  --interval-padding 100 \
  --annovar-bin  /path/to/annovar \
  --annovar-db   /path/to/annovar/humandb \
  --threads      32 \
  --fastp-threads 8 \
  --annovar-threads 24 \
  --max-workers  2 \
  --java-opts    "-Xmx80g"
```

### V2 modes

```bash
# Cohort joint genotyping — one shared VCF/annotation instead of per-sample files
exomeflow run --input-dir fastq/ --output results/ --joint-genotyping --intervals targets.bed

# Somatic tumor-only calling
exomeflow run --input-dir fastq/ --output results/ --mode somatic --germline-resource af-only-gnomad.vcf.gz

# Read-depth CNV alongside germline calling
exomeflow run --input-dir fastq/ --output results/ --cnv --intervals targets.bed

# GRCh37/hg19
exomeflow run --input-dir fastq/ --output results/ --genome-build GRCh37
```

Every sample gets its own separate annotated output file by default, no matter how
many samples are in `--input-dir` — cohort output only happens with `--joint-genotyping`,
and that's an explicit opt-in, never an automatic side effect of batching samples.

### Check version

```bash
exomeflow --version
```

### Get help

```bash
exomeflow run --help
```

---

## 5. All Command-Line Options

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `--input-dir` | Yes | — | Directory containing paired FASTQ files |
| `--output` | Yes | `results/` | Root output directory |
| `--reference` | No | auto-resolved | BWA-indexed reference genome FASTA |
| `--dbsnp` | No | auto-resolved | dbSNP VCF (bgzipped + tabix-indexed) |
| `--mills` | No | auto-resolved | Mills and 1000G gold standard indels VCF |
| `--known-indels` | No | auto-resolved | Known indels VCF |
| `--intervals` | No | — | Exome capture BED. Required for `--joint-genotyping`/`--cnv` |
| `--interval-padding` | No | `100` | Base-pair padding around each target interval |
| `--annovar-bin` | No | auto-resolved | ANNOVAR installation directory |
| `--annovar-db` | No | auto-resolved | ANNOVAR humandb directory |
| `--mode` | No | `germline` | `germline` (HaplotypeCaller) or `somatic` (tumor-only Mutect2) |
| `--genome-build` | No | `hg38` | `hg38` or `GRCh37` |
| `--joint-genotyping` | No | off | Cohort mode — one shared VCF/annotation, not per-sample |
| `--cnv` | No | off | Also call read-depth CNVs per sample |
| `--germline-resource` | No | — | gnomAD AF-only VCF for `--mode somatic` |
| `--threads` | No | `24` | Threads for BWA MEM and HaplotypeCaller |
| `--fastp-threads` | No | `8` | Threads for fastp |
| `--annovar-threads` | No | `24` | Threads for ANNOVAR |
| `--max-workers` | No | `1` | Number of samples to process in parallel |
| `--java-opts` | No | `-Xmx80g` | JVM memory options passed to GATK |
| `--version` | No | — | Print version and exit |

First run auto-resolves the six "auto-resolved" paths above interactively (or
downloads them) and saves the result to `~/.exomeflow/config.json` — pass any of them
explicitly to override the saved value for a single run.

---

## 6. Understanding the Output

After a successful run, the `results/` directory contains:

```
results/
│
├── filtered_fastp/                         ← Step 1: fastp QC output
│   ├── sample1_1_filtered.fastq.gz
│   ├── sample1_2_filtered.fastq.gz
│   ├── sample1_fastp.html                  ← Open in browser for QC report
│   └── sample1_fastp.json
│
├── Mapsam/                                 ← Steps 2–7: Alignment + BQSR
│   ├── sample1_recalibrated.bam            ← Final analysis-ready BAM
│   ├── sample1_recalibrated.bam.bai        ← BAM index
│   └── sample1_flagstat.txt                ← Alignment statistics
│
├── VCF/                                    ← Variant calling (per-sample by default)
│   ├── sample1.vcf                         ← Raw HaplotypeCaller output
│   ├── sample1_PASS.vcf                    ← PASS-only filtered variants
│   ├── sample1.annovar.hg38_multianno.vcf  ← Annotated VCF
│   ├── sample1.annovar.hg38_multianno.txt  ← Annotated tab-delimited table
│   ├── sample1.annovar.hpo.txt             ← + HPO terms + ACMG classification
│   └── cohort/                             ← Only with --joint-genotyping
│       ├── cohort.vcf.gz
│       ├── cohort_PASS.vcf
│       └── cohort.annovar.hg38_multianno.{vcf,txt}
│
├── CNV/                                    ← Only with --cnv
│   └── sample1_denoised_cr.tsv + plot
│
├── multiqc/
│   └── exomeflow_report.html               ← Cohort-wide QC rollup
│
├── logs/
│   ├── analysis_20250101_120000.log        ← Full pipeline log
│   ├── errors_20250101_120000.log          ← Errors only
│   └── sample1_20250101_120000.log         ← Per-sample log
│
└── .checkpoints/                           ← Resume state (do not delete)
    ├── sample1.fastp.done
    ├── sample1.bwa.done
    └── ...
```

### Key output files explained

| File | Use |
|------|-----|
| `*_recalibrated.bam` | Load into **IGV** to visually inspect variants |
| `*_PASS.vcf` | Clean variant list — use for downstream analysis |
| `*.hg38_multianno.txt` | Open in Excel / R for variant interpretation |
| `*.hg38_multianno.vcf` | Annotated VCF — submit to clinical databases |
| `*_fastp.html` | QC report — check read quality before trusting results |
| `*_flagstat.txt` | Check alignment rate — should be > 95% |

---

## 7. Resuming an Interrupted Run

If the pipeline is interrupted (power cut, time limit, crash), simply
**run the exact same command again**. ExomeFlow checks `.checkpoints/`
and skips every step that already completed.

```bash
# First run — interrupted at BQSR
exomeflow run --input-dir fastq/ --output results/ ...

# Resume — automatically skips fastp, BWA, sort, flagstat, markdup, index
exomeflow run --input-dir fastq/ --output results/ ...
```

> **Important:** Do not delete the `results/.checkpoints/` directory
> between runs if you want to resume.

To force a full re-run from scratch:

```bash
rm -rf results/.checkpoints/
exomeflow run --input-dir fastq/ --output results/ ...
```

---

## 8. Running Multiple Samples in Parallel

By default ExomeFlow processes one sample at a time (`--max-workers 1`).  
To process multiple samples simultaneously:

```bash
exomeflow run \
  --input-dir fastq/ \
  --output    results/ \
  --max-workers 4 \
  --threads 8 \
  ...
```

### How to choose `--max-workers` and `--threads`

```
Total CPU cores = max_workers × threads

Example — 48-core server:
  --max-workers 2 --threads 24   → 2 samples × 24 threads = 48 cores used
  --max-workers 4 --threads 12   → 4 samples × 12 threads = 48 cores used
  --max-workers 1 --threads 48   → 1 sample  × 48 threads = 48 cores used
```

> **Rule of thumb:** For a server with 48+ cores and large RAM (> 64 GB),
> `--max-workers 2 --threads 24` is a good starting point.

---

## 9. Common Errors and Fixes

### Requirements check fails at startup

```
EnvironmentError: Requirements check failed — fix the following issues:
  • gatk not found on PATH → add to PATH: export PATH=/path/to/gatk-4.6.x.x:$PATH
```

**Fix:** Add the missing tool to your PATH and re-run.

---

### No FASTQ files found

```
FileNotFoundError: No paired FASTQ files (matching *_1.fastq.gz or *_R1.fastq.gz) found in fastq/
```

**Fix:** Check that your files follow one of the two accepted naming conventions:
`<sample_id>_1.fastq.gz` / `<sample_id>_2.fastq.gz`, or
`<sample_id>_R1.fastq.gz` / `<sample_id>_R2.fastq.gz`

---

### GATK runs out of memory

```
java.lang.OutOfMemoryError: Java heap space
```

**Fix:** Increase the JVM heap size:

```bash
exomeflow run ... --java-opts "-Xmx120g"
```

---

### HaplotypeCaller is very slow

**Cause:** You did not provide an exome intervals BED file — the pipeline
falls back to whole-genome mode.

**Fix:** Add `--intervals` pointing to your capture kit BED:

```bash
exomeflow run ... --intervals /path/to/exome_targets.bed
```

---

### ANNOVAR annotation fails

```
PipelineStepError: table_annovar.pl failed
```

**Fixes to check:**
1. Verify `--annovar-bin` points to the directory containing `table_annovar.pl`
2. Verify `--annovar-db` contains the hg38 database files
3. Make sure the required databases are downloaded:

```bash
perl /path/to/annovar/annotate_variation.pl \
  -buildver hg38 \
  -downdb -webfrom annovar refGene /path/to/annovar/humandb/

perl /path/to/annovar/annotate_variation.pl \
  -buildver hg38 \
  -downdb -webfrom annovar clinvar_20240416 /path/to/annovar/humandb/
```

---

### A sample failed but others completed

ExomeFlow logs the failure and **continues processing other samples**.  
Check the per-sample log for details:

```bash
cat results/logs/sample1_*.log | grep ERROR
```

The failed sample can be re-run after fixing the issue — checkpointing
ensures all completed steps are skipped.

---

## Quick Reference Card

```bash
# Install
pip install exomeflow

# Check requirements
python check_requirements.py

# Run (single sample, 24 threads)
exomeflow run \
  --input-dir fastq/ --output results/ \
  --reference refs/hg38.fa \
  --dbsnp refs/dbsnp.vcf.gz \
  --mills refs/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known-indels refs/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --intervals refs/exome_targets.bed \
  --annovar-bin /path/to/annovar \
  --annovar-db /path/to/annovar/humandb \
  --threads 24

# Run (cohort, 2 samples in parallel)
exomeflow run ... --max-workers 2 --threads 24

# Resume after interruption
exomeflow run ...   # (same command — checkpoints skip completed steps)

# Get help
exomeflow run --help
exomeflow --version
```

---

*ExomeFlow v2.1.2 — Robin Kumar, AIIMS New Delhi, 2026*
