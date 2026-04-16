# ============================================================
# ExomeFlow — Production WES Analysis Pipeline
# Docker image
#
# Base: broadinstitute/gatk:4.6.2.0 (official GATK image)
#       Includes: GATK 4.6.2.0, Java 17, Miniconda3, Ubuntu
#
# Build:
#   docker build -t itsrobintomar/exomeflow:1.0.6 .
#   docker tag itsrobintomar/exomeflow:1.0.6 itsrobintomar/exomeflow:latest
#
# Push:
#   docker push itsrobintomar/exomeflow:1.0.6
#   docker push itsrobintomar/exomeflow:latest
#
# Run:
#   docker run --rm -it \
#     -v /path/to/fastq:/data/fastq \
#     -v /path/to/refs:/refs \
#     -v /path/to/annovar:/annovar \
#     -v /path/to/results:/data/results \
#     itsrobintomar/exomeflow:latest run \
#       --input-dir    /data/fastq \
#       --output       /data/results \
#       --reference    /refs/hg38.fa \
#       --dbsnp        /refs/dbsnp.vcf.gz \
#       --mills        /refs/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
#       --known-indels /refs/Homo_sapiens_assembly38.known_indels.vcf.gz \
#       --annovar-bin  /annovar \
#       --annovar-db   /annovar/humandb \
#       --threads      24
#
# NOTE: ANNOVAR requires registration at annovar.openbioinformatics.org
#       Mount your local ANNOVAR installation: -v /path/to/annovar:/annovar
# ============================================================

# Official Broad Institute GATK image — GATK 4.6.2.0 + Java 17 + Miniconda pre-installed
FROM broadinstitute/gatk:4.6.2.0

LABEL maintainer="Robin Tomar <itsrobintomar@gmail.com>"
LABEL org.opencontainers.image.title="ExomeFlow"
LABEL org.opencontainers.image.description="Production-quality Whole Exome Sequencing analysis pipeline"
LABEL org.opencontainers.image.version="1.0.6"
LABEL org.opencontainers.image.url="https://github.com/imrobintomar/exomeflow"
LABEL org.opencontainers.image.licenses="MIT"

# ── System packages ──────────────────────────────────────────────────────────
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        perl \
        libncurses5-dev \
        libbz2-dev \
        liblzma-dev \
        zlib1g-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ── Conda channels ───────────────────────────────────────────────────────────
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority flexible

# ── Bioinformatics tools (GATK already installed in base image) ───────────────
RUN conda install -y \
        "bwa>=0.7.17" \
        "samtools>=1.13" \
        "fastp>=0.20.1" \
        perl \
    && conda clean -afy

# ── Java options ─────────────────────────────────────────────────────────────
ENV JAVA_OPTS="-Xmx80g"

# ── ExomeFlow Python package ─────────────────────────────────────────────────
RUN pip install --no-cache-dir exomeflow==1.0.6

# ── Verify all tools ─────────────────────────────────────────────────────────
RUN echo "=== Verifying tools ===" && \
    gatk --version && \
    bwa 2>&1 | grep -i "version\|Program" | head -2 && \
    samtools --version | head -1 && \
    fastp --version 2>&1 | head -1 && \
    perl --version | head -2 && \
    exomeflow --version && \
    echo "=== All tools verified ==="

# ── Mount points ─────────────────────────────────────────────────────────────
# /data/fastq    — input FASTQ files
# /data/results  — pipeline output
# /refs          — reference genome + VCF files (hg38.fa, dbsnp.vcf.gz, etc.)
# /annovar       — ANNOVAR installation (must be mounted by user)
VOLUME ["/data/fastq", "/data/results", "/refs", "/annovar"]

WORKDIR /data

ENTRYPOINT ["exomeflow"]
CMD ["--help"]
