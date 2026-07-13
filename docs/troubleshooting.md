# Troubleshooting

## No FASTQ files found

```
No paired FASTQ files (matching *_1.fastq.gz) found in fastq/
```

This is checked before the setup wizard runs, so it fails in under a second
rather than after a multi-hour first-run download. Check that your files
follow the naming convention: `<sample_id>_1.fastq.gz` and
`<sample_id>_2.fastq.gz`, both in `--input-dir`.

## `--mode` rejected

```
--mode must be 'germline' or 'somatic', got '<value>'.
```

Only those two values are valid for `--mode`.

## `--joint-genotyping` / `--cnv` without `--intervals`

```
--joint-genotyping and --cnv both require --intervals
(a bounded region is needed before either can run).
```

Both cohort joint genotyping and CNV calling need a bounded target region.
Pass `--intervals /path/to/capture_kit.bed`.

## `--annovar-protocols` / `--annovar-operations` mismatch

```
--annovar-protocols and --annovar-operations must be given together
(ANNOVAR pairs each protocol with an operation).
```

or

```
--annovar-protocols and --annovar-operations have a different number
of comma-separated entries.
```

These two flags exist to let you work around a pre-existing, differently
versioned ANNOVAR install (e.g. a newer/older ClinVar build) without editing
source — pass both together with matching entry counts, or neither.

## GATK runs out of memory

```
java.lang.OutOfMemoryError: Java heap space
```

`--java-opts` is auto-sized from available RAM (60% of available, capped
4–80g) if you don't pass it explicitly. On a memory-constrained machine,
increase it manually:

```bash
exomeflow run ... --java-opts "-Xmx120g"
```

## HaplotypeCaller is very slow

**Cause:** no `--intervals` was given, so the pipeline calls across the whole
genome instead of just the exome capture region.

**Fix:**

```bash
exomeflow run ... --intervals /path/to/exome_targets.bed
```

## ANNOVAR annotation fails

```
PipelineStepError: table_annovar.pl failed
```

Check:

1. `--annovar-bin` points to the directory containing `table_annovar.pl`.
2. `--annovar-db` contains the database files for your `--genome-build`.
3. The required databases are actually downloaded — ExomeFlow checks
   per-database-file completeness on every run and will re-fetch anything
   missing, but a network failure mid-download can leave a partial database
   that needs re-running `exomeflow setup` to repair.

## ANNOVAR / annotation skipped with a warning on an empty VCF

If a sample's `_PASS.vcf` has zero variants (e.g. a catastrophic sequencing
failure, or an intervals BED with near-zero coverage overlap), annotation is
skipped cleanly with a warning instead of crashing — ANNOVAR itself cannot
process a header-only VCF. This is expected behavior for that specific
input, not a bug.

## A sample failed but others completed

ExomeFlow logs the failure and continues processing other samples in the
same run. Check the per-sample log for details:

```bash
grep ERROR results/logs/sample1_*.log
```

Fix the underlying cause, then re-run the exact same command — checkpointing
in `results/.checkpoints/` skips every step that already completed and
retries the failed sample and everything after it.

## Switching `--genome-build` on an existing output directory

Changing `--genome-build` between runs against the same `--output` is
supported — ExomeFlow detects the build change against the previously saved
config and re-resolves reference/ANNOVAR paths for the new build rather than
silently reusing stale hg38 files against a GRCh37 run (or vice versa).
