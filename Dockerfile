# ============================================================
# ExomeFlow — Production WES Analysis Pipeline
# Public image — safe to push to Docker Hub
#
# Base  : broadinstitute/gatk:4.6.2.0
#         Includes GATK 4.6.2.0 · Java 17 · Miniconda3 · Ubuntu
#
# ANNOVAR: cannot be bundled (license). Mount your local copy at runtime.
#
# ── Build & push ────────────────────────────────────────────
#   docker build -t itsrobintomar/exomeflow:2.2.11 .
#   docker tag  itsrobintomar/exomeflow:2.2.11 itsrobintomar/exomeflow:latest
#   docker push itsrobintomar/exomeflow:2.2.11
#   docker push itsrobintomar/exomeflow:latest
#
# ── Run ─────────────────────────────────────────────────────
#   docker run --rm -it \
#     -v /path/to/fastq:/data/fastq \
#     -v /path/to/refs:/refs \
#     -v /path/to/vcf:/vcf \
#     -v /path/to/annovar:/annovar \
#     -v /path/to/results:/data/results \
#     itsrobintomar/exomeflow:latest run \
#       --input-dir    /data/fastq \
#       --output       /data/results \
#       --reference    /refs/Homo_sapiens_assembly38.fasta \
#       --dbsnp        /vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
#       --mills        /vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
#       --known-indels /vcf/Homo_sapiens_assembly38.known_indels.vcf.gz \
#       --annovar-bin  /annovar \
#       --annovar-db   /annovar/humandb \
#       --threads      24
# ============================================================

FROM broadinstitute/gatk:4.6.2.0

LABEL maintainer="Robin Kumar <itsrobintomar@gmail.com>"
LABEL org.opencontainers.image.title="ExomeFlow"
LABEL org.opencontainers.image.description="Production-quality Whole Exome Sequencing analysis pipeline"
LABEL org.opencontainers.image.version="2.2.11"
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

# ── Bioinformatics tools ─────────────────────────────────────────────────────
# GATK already in base image — only add bwa, samtools, fastp, perl
RUN conda install -y \
        "bwa>=0.7.17" \
        "samtools>=1.13" \
        "fastp>=0.20.1" \
        perl \
    && conda clean -afy

# ── JVM options ──────────────────────────────────────────────────────────────
ENV JAVA_OPTS="-Xmx80g"

# ── ExomeFlow ────────────────────────────────────────────────────────────────
RUN pip install --no-cache-dir exomeflow==2.2.11

# ── Verify ───────────────────────────────────────────────────────────────────
RUN gatk --version && \
    bwa 2>&1 | head -3 && \
    samtools --version | head -1 && \
    fastp --version 2>&1 | head -1 && \
    perl --version | head -2 && \
    exomeflow --version

# ── Volumes ──────────────────────────────────────────────────────────────────
# /data/fastq    — input FASTQ files  (*_1.fastq.gz / *_2.fastq.gz)
# /data/results  — pipeline output
# /refs          — reference FASTA + BWA index (Homo_sapiens_assembly38.fasta)
# /vcf           — dbsnp138, Mills indels, known indels VCFs
# /annovar       — ANNOVAR scripts (table_annovar.pl etc.) — USER MUST MOUNT
# /annovar/humandb — ANNOVAR databases — USER MUST MOUNT
VOLUME ["/data/fastq", "/data/results", "/refs", "/vcf", "/annovar"]

WORKDIR /data

ENTRYPOINT ["exomeflow"]
CMD ["--help"]
