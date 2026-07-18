# Modes & Options

ExomeFlow's default behavior (germline calling, hg38, per-sample output) never
changes unless you opt in to one of the flags below. Every capability here
plugs into the same pluggable step registry described in
[Architecture](architecture.md) — turning a flag on/off adds or removes steps
from the pipeline, it doesn't fork into a separate code path.

## Cohort joint genotyping

```bash
exomeflow run --input-dir fastq/ --output results/ \
  --joint-genotyping --intervals targets.bed
```

Runs HaplotypeCaller in `-ERC GVCF` mode per sample, then GenomicsDBImport →
GenotypeGVCFs across all samples, producing one shared cohort VCF and one
shared annotated file instead of per-sample output:

```
VCF/cohort/cohort.vcf.gz
VCF/cohort/cohort_PASS.vcf
VCF/cohort/cohort.annovar.hg38_multianno.{vcf,txt}
```

Requires `--intervals` — joint genotyping without a target region set is not
supported. This is strictly opt-in: running any number of samples without this
flag still produces one annotated file per sample.

## Somatic mode (tumor-only)

```bash
exomeflow run --input-dir fastq/ --output results/ \
  --mode somatic \
  --germline-resource af-only-gnomad.vcf.gz \
  --panel-of-normals pon.vcf.gz
```

Replaces HaplotypeCaller with Mutect2 in tumor-only mode, and the germline
hard-filter chain with GATK's contamination-aware filtering chain
(GetPileupSummaries → CalculateContamination → FilterMutectCalls). Tumor-normal
paired calling is not yet supported — `--germline-resource` (a gnomAD AF-only
VCF) and `--panel-of-normals` (a pre-built PoN VCF) are both optional but
strongly recommended: this is GATK's own documented alternative to matched
tumor-normal calling, not a lesser substitute. Without them, tumor-only calls
carry more false positives and aren't filtered against recurrent sequencing
artifacts a PoN would catch. ExomeFlow consumes an existing PoN VCF; building
one from a set of normal samples isn't automated yet — see GATK's
`CreateSomaticPanelOfNormals` workflow.

## Read-depth CNV calling

```bash
exomeflow run --input-dir fastq/ --output results/ \
  --cnv --intervals targets.bed
```

Adds GATK's read-depth CNV chain (CollectReadCounts → DenoiseReadCounts →
PlotDenoisedCopyRatios) per sample, alongside — not instead of — SNP/INDEL
calling. Requires `--intervals`. This runs in single-sample mode without a
panel of normals, which is less accurate than a proper cohort PoN — treat
CNV calls as a starting point for review, not a final answer.

```
CNV/sample1_denoised_cr.tsv
CNV/sample1_denoised_cr.png
```

## GRCh37 / hg19

```bash
exomeflow run --input-dir fastq/ --output results/ --genome-build GRCh37
```

Reference bundle, ANNOVAR buildver, and database downloads all follow the
selected build automatically. One caveat worth knowing: GRCh37 annotation uses
`gnomad211_exome` / `gnomad211_genome` instead of hg38's
`gnomad41_exome` / `gnomad41_genome` — gnomAD v4.1 was never released for
hg19/GRCh37, so this substitution is intentional, not a downgrade you can fix
by asking for a newer database.

## HPO terms + ACMG classification

On by default, no flag required. After ANNOVAR annotation, every sample's
`*.annovar.hg38_multianno.txt` is joined with HPO gene-to-phenotype terms and
run through InterVar for ACMG/AMP classification, written to
`*.annovar.hpo.txt` alongside the raw ANNOVAR output. If HPO data or InterVar
can't be auto-provisioned (e.g. no network on first run), this step degrades
to a skip-with-warning rather than failing the pipeline — the raw annotated
output is still produced either way.

## MultiQC rollup

Also on by default. After all samples finish, ExomeFlow generates
`multiqc/exomeflow_report.html` — a single HTML report rolling up fastp, BWA,
MarkDuplicates, and BQSR metrics across every sample in the run.
