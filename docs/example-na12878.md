# Real-World Example — NA12878

This walkthrough runs ExomeFlow 2.1.1 end-to-end against
[NA12878/HG001](https://www.coriell.org/0/Sections/Search/Sample_Detail.aspx?Ref=GM12878),
the Genome in a Bottle consortium's most widely used reference sample, using a
real commercial exome capture kit's target region file (Agilent SureSelect
V8 CES) — not a truth-set confidence BED, so the numbers below reflect what a
typical clinical exome run actually looks like, not an artificially narrowed
benchmarking region.

## The command

```bash
exomeflow run \
  --input-dir  benchmarks/Samples \
  --output     benchmarks/results \
  --intervals  hg38_agilentV8_CES.bed \
  --genome-build hg38
```

On a machine that already has the setup wizard's dependencies resolved, this
is genuinely the whole command — no reference paths, no ANNOVAR database
paths, nothing else to configure.

!!! note "Numbers on this page are being refreshed"
    A full re-run against NA12878 on the current 2.1.1 release is in progress
    (using a real exome-capture-kit BED, not a truth-set confidence region) to
    replace the walkthrough below with real variant counts, runtime, and
    accuracy figures. This page will be updated in place once it completes —
    check back soon, or see the [changelog](changelog.md) for what's new in
    2.1.1.

## What to expect while it runs

```
[NA12878] Running fastp ...
[NA12878] Running BWA MEM ...
[NA12878] Sorting BAM ...
[NA12878] Running flagstat ...
[NA12878] Marking duplicates ...
[NA12878] Indexing BAM ...
[NA12878] Running BQSR ...
[NA12878] Running HaplotypeCaller ...
[NA12878] Filtering variants ...
[NA12878] Running ANNOVAR annotation ...
[NA12878] Joining HPO terms ...
[NA12878] Running ACMG classification (InterVar) ...
Generating MultiQC report ...
```

Alignment and variant calling dominate the runtime; ANNOVAR annotation and the
HPO/ACMG enrichment step are comparatively fast once the databases are already
downloaded.

## Output

```
results/VCF/NA12878.vcf                          ← raw HaplotypeCaller calls
results/VCF/NA12878_PASS.vcf                      ← PASS-only variants
results/VCF/NA12878.annovar.hg38_multianno.txt    ← annotated table
results/VCF/NA12878.annovar.hpo.txt               ← + HPO terms + ACMG class
results/Mapsam/NA12878_flagstat.txt               ← alignment stats
results/multiqc/exomeflow_report.html             ← QC rollup
```

See [Output Files](output-files.md) for what every file in this layout
contains.
