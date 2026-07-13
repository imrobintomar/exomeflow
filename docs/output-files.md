# Output Files

After a successful run, `--output` contains:

```
results/
│
├── filtered_fastp/                         ← fastp QC output
│   ├── sample1_1_filtered.fastq.gz
│   ├── sample1_2_filtered.fastq.gz
│   ├── sample1_fastp.html                  ← Open in browser for QC report
│   └── sample1_fastp.json
│
├── Mapsam/                                 ← Alignment + BQSR
│   ├── sample1_recalibrated.bam            ← Final analysis-ready BAM
│   ├── sample1_recalibrated.bam.bai
│   └── sample1_flagstat.txt                ← Alignment statistics
│
├── VCF/                                    ← Variant calling (per-sample by default)
│   ├── sample1.vcf                         ← Raw HaplotypeCaller/Mutect2 output
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
│   └── sample1_denoised_cr.tsv + .png
│
├── multiqc/
│   └── exomeflow_report.html               ← Cohort-wide QC rollup
│
├── logs/
│   ├── analysis_<timestamp>.log            ← Full pipeline log
│   ├── errors_<timestamp>.log              ← Errors only
│   └── sample1_<timestamp>.log             ← Per-sample log
│
└── .checkpoints/                           ← Resume state — do not delete
    ├── sample1.fastp.done
    ├── sample1.bwa.done
    └── ...
```

## Key files at a glance

| File | Use |
|---|---|
| `*_recalibrated.bam` | Load into IGV to visually inspect variants. |
| `*_PASS.vcf` | Clean variant list — use for downstream analysis. |
| `*.hg38_multianno.txt` | Open in Excel / R for variant interpretation. |
| `*.annovar.hpo.txt` | Same table, joined with HPO terms + ACMG classification. |
| `*.hg38_multianno.vcf` | Annotated VCF — submit to clinical databases. |
| `*_fastp.html` | QC report — check read quality before trusting downstream results. |
| `*_flagstat.txt` | Check alignment rate — should be > 95% for a healthy exome run. |
| `exomeflow_report.html` | One-page rollup across every sample in the run. |

See [Real-World Example — NA12878](example-na12878.md) for what these files
actually contain on a real run.
