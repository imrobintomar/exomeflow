# CLI Reference

This page mirrors `exomeflow run --help` on the currently published version. Run
it yourself any time to get the exact reference for your installed version.

## `exomeflow run`

```
exomeflow run --input-dir fastq/ --output results/
```

### Default workflow (germline, single-sample)

```
FASTQ → fastp → BWA MEM → BAM processing → BQSR → HaplotypeCaller
     → hard-filter → ANNOVAR → HPO terms → ACMG classification
...then once, across all samples: MultiQC rollup.
```

Every sample gets its own separate annotated output file by default, no matter
how many samples are in `--input-dir` — cohort output only happens with
`--joint-genotyping`, and that's an explicit opt-in, never an automatic side
effect of batching samples.

### Options

| Option | Default | Description |
|---|---|---|
| `--input-dir`, `-i` | *(required)* | Directory containing paired FASTQ files (`*_1.fastq.gz` / `*_2.fastq.gz`). |
| `--output`, `-o` | `results` | Root output directory (created if absent). |
| `--reference`, `-r` | auto-resolved | Reference genome FASTA (hg38 or GRCh37). |
| `--dbsnp` | auto-resolved | dbSNP VCF. |
| `--mills` | auto-resolved | Mills and 1000G gold standard indels VCF. |
| `--known-indels` | auto-resolved | Known indels VCF. |
| `--annovar-bin` | auto-resolved | ANNOVAR installation directory. |
| `--annovar-db` | auto-resolved | ANNOVAR humandb directory. |
| `--annovar-protocols` | build default | Override the ANNOVAR `--protocol` list (comma-separated). |
| `--annovar-operations` | build default | Override the ANNOVAR `--operation` list — must match `--annovar-protocols` in count. |
| `--intervals` | — | Exome capture BED file. Omit for whole-genome mode. |
| `--interval-padding` | `100` | Base-pair padding around each target interval. |
| `--threads`, `-t` | auto-detected | Threads for BWA MEM and HaplotypeCaller. |
| `--fastp-threads` | `8` | Threads for fastp. |
| `--annovar-threads` | `24` | Threads for ANNOVAR. |
| `--max-workers` | `1` | Number of samples to process in parallel. |
| `--java-opts` | auto-sized | JVM options passed via `JAVA_OPTS` (60% of available RAM, 4–80g). |
| `--mode` | `germline` | `germline` (HaplotypeCaller) or `somatic` (tumor-only Mutect2). |
| `--genome-build` | `hg38` | `hg38` or `GRCh37`. |
| `--joint-genotyping` | off | Cohort mode: one shared VCF instead of per-sample files. |
| `--cnv` | off | Also call read-depth CNVs per sample (requires `--intervals`). |
| `--germline-resource` | — | gnomAD AF-only VCF for `--mode somatic` (optional but recommended). |
| `--panel-of-normals` | — | Pre-built Panel of Normals VCF for `--mode somatic` (optional but recommended). |
| `--yes`, `-y` | off | Non-interactive: auto-confirm every setup prompt. |

Any of the auto-resolved paths can be passed explicitly to override the value
saved in `~/.exomeflow/config.json` for a single run.

## `exomeflow setup`

Re-runs the setup wizard directly — change a reference path, point at a
different ANNOVAR database directory, or repair a partially-provisioned config
without running the full pipeline.

```bash
exomeflow setup
```

## `exomeflow --version`

```bash
exomeflow --version
```

## Resuming an interrupted run

If the pipeline is interrupted (power cut, time limit, crash), run the exact
same command again. ExomeFlow checks `results/.checkpoints/` and skips every
step that already completed for every sample — don't delete that directory
between runs if you want to resume.

```bash
rm -rf results/.checkpoints/   # force a full re-run from scratch instead
```

## Running multiple samples in parallel

```bash
exomeflow run --input-dir fastq/ --output results/ --max-workers 4 --threads 8
```

```
Total CPU cores = max_workers × threads

Example — 48-core server:
  --max-workers 2 --threads 24   → 2 samples × 24 threads = 48 cores used
  --max-workers 4 --threads 12   → 4 samples × 12 threads = 48 cores used
  --max-workers 1 --threads 48   → 1 sample  × 48 threads = 48 cores used
```
