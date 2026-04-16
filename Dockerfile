# ============================================================
# ExomeFlow — Production WES Analysis Pipeline
# Docker image
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
#       --annovar-db   /annovar/humandb
#
# NOTE: ANNOVAR requires registration at annovar.openbioinformatics.org
#       Mount your local ANNOVAR installation using -v /path/to/annovar:/annovar
# ============================================================

FROM continuumio/miniconda3:24.1.2-0

LABEL maintainer="Robin Tomar <itsrobintomar@gmail.com>"
LABEL org.opencontainers.image.title="ExomeFlow"
LABEL org.opencontainers.image.description="Production-quality Whole Exome Sequencing analysis pipeline"
LABEL org.opencontainers.image.version="1.0.6"
LABEL org.opencontainers.image.url="https://github.com/itsrobintomar/exomeflow"
LABEL org.opencontainers.image.licenses="MIT"

# ── System packages ──────────────────────────────────────────────────────────
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        wget curl git perl \
        libncurses5-dev libbz2-dev liblzma-dev \
        zlib1g-dev libssl-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ── Conda channels ───────────────────────────────────────────────────────────
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority flexible

# ── Bioinformatics tools ─────────────────────────────────────────────────────
RUN conda install -y \
        "bwa>=0.7.17" \
        "samtools>=1.13" \
        "fastp>=0.20.1" \
        "gatk4>=4.4.0" \
        perl \
        "openjdk>=17" \
    && conda clean -afy

# ── Java options ─────────────────────────────────────────────────────────────
ENV JAVA_OPTS="-Xmx80g"

# ── ExomeFlow Python package ─────────────────────────────────────────────────
RUN pip install --no-cache-dir exomeflow==1.0.6

# ── Verify installation ──────────────────────────────────────────────────────
RUN exomeflow --version && \
    bwa 2>&1 | head -1 && \
    samtools --version | head -1 && \
    fastp --version 2>&1 | head -1 && \
    gatk --version

# ── Mount points ─────────────────────────────────────────────────────────────
# /data/fastq    — input FASTQ files
# /data/results  — pipeline output
# /refs          — reference genome + VCF files
# /annovar       — ANNOVAR installation (must be mounted by user)
VOLUME ["/data/fastq", "/data/results", "/refs", "/annovar"]

WORKDIR /data

ENTRYPOINT ["exomeflow"]
CMD ["--help"]
