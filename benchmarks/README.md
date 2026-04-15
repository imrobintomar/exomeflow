# ExomeFlow — Benchmarking Suite

Four standalone scripts. Run them after the pipeline has processed your samples.

```
benchmarks/
├── 01_performance.py      Runtime per step (parses log files)
├── 02_variant_stats.py    Ts/Tv, Het/Hom, dbSNP concordance (bcftools stats)
├── 03_annovar_stats.py    Functional annotation stats (ANNOVAR multianno.txt)
└── 04_accuracy.py         Precision / Recall / F1 vs GIAB truth set (hap.py)
```

---

## Install Python dependencies

```bash
pip install pandas matplotlib seaborn numpy
```

---

## Script 1 — Performance (runtime per step)

Parses ExomeFlow log files. No extra tools needed.

```bash
python benchmarks/01_performance.py \
    --log-dir results/logs \
    --out-dir benchmarks/results
```

| Output | Description |
|--------|-------------|
| `runtime_table.csv` | Step × Sample runtime matrix (minutes) |
| `runtime_heatmap.png` | Colour heatmap of runtimes |
| `runtime_bars.png` | Stacked bar chart per sample |

---

## Script 2 — Variant quality statistics

Runs `bcftools stats` on every `*_PASS.vcf`. Works on all samples without a truth set.

```bash
python benchmarks/02_variant_stats.py \
    --vcf-dir results/VCF \
    --dbsnp   /data/refs/dbsnp.vcf.gz \
    --out-dir benchmarks/results
```

| Output | Description |
|--------|-------------|
| `variant_stats.csv` | SNPs, INDELs, Ts/Tv, Het/Hom, dbSNP % per sample |
| `variant_counts.png` | SNP/INDEL bar chart |
| `tstv_hethom.png` | QC plots with expected ranges shaded |

**Expected ranges for a good WES run**

| Metric | Expected |
|--------|----------|
| Ts/Tv ratio | 2.0 – 3.3 |
| Het/Hom ratio | 1.5 – 2.5 |
| dbSNP concordance | > 95% |

---

## Script 3 — ANNOVAR functional annotation statistics

Parses `*.annovar.hg38_multianno.txt` files. No extra tools needed.

```bash
python benchmarks/03_annovar_stats.py \
    --vcf-dir results/VCF \
    --out-dir benchmarks/results
```

| Output | Description |
|--------|-------------|
| `annovar_stats.csv` | Per-sample counts: exonic, splicing, pathogenic, novel |
| `func_distribution.png` | Functional consequence stacked bar chart |
| `clinvar_distribution.png` | ClinVar classification per sample |
| `gnomad_af.png` | gnomAD allele frequency spectrum |
| `novel_vs_known.png` | Novel (not in dbSNP) vs known variants |

---

## Script 4 — Accuracy benchmarking (GIAB truth set)

Compares your calls against NA12878 / HG001 GIAB truth using hap.py.

### Step 1 — Install hap.py
```bash
conda install -c bioconda hap.py
```

### Step 2 — Download GIAB NA12878 truth files (~1 GB)
```bash
# Create a directory
mkdir -p /data/refs/giab && cd /data/refs/giab

# Truth VCF + index
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

# High-confidence BED regions
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
```

### Step 3 — Download NA12878 FASTQ (if not already in your dataset)
```bash
# From SRA (~50 GB, 30x WGS)
fastq-dump --split-files --gzip SRR622461
# Rename to match ExomeFlow convention
mv SRR622461_1.fastq.gz NA12878_1.fastq.gz
mv SRR622461_2.fastq.gz NA12878_2.fastq.gz
```

### Step 4 — Run accuracy benchmark
```bash
python benchmarks/04_accuracy.py \
    --query-vcf  results/VCF/NA12878_PASS.vcf \
    --truth-vcf  /data/refs/giab/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    --truth-bed  /data/refs/giab/HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
    --reference  /data/refs/hg38.fa \
    --sample-name NA12878
```

| Output | Description |
|--------|-------------|
| `accuracy_table.csv` | Precision / Recall / F1 for SNP and INDEL |
| `accuracy_snp_indel.png` | Bar chart of all three metrics |
| `precision_recall.png` | Precision vs Recall scatter plot |

**Expected values for a well-run GATK pipeline vs GIAB**

| Variant type | Metric | Expected |
|---|---|---|
| SNP | Precision | > 0.99 |
| SNP | Recall | > 0.99 |
| SNP | F1 | > 0.99 |
| INDEL | Precision | > 0.95 |
| INDEL | Recall | > 0.95 |
| INDEL | F1 | > 0.95 |
