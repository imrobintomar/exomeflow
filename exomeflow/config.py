"""
Pipeline configuration dataclass.

All paths and tuning parameters are held here.  CLI arguments are merged
into a Config instance before the pipeline starts so every module receives
a single, consistent source of truth.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path


# ANNOVAR annotation databases — mirrors the Bash script's constants
ANNOVAR_PROTOCOLS: str = (
    "refGene,dbnsfp47a,clinvar_20240416,gnomad41_exome,"
    "gnomad41_genome,avsnp150,cosmic84_coding,exac03"
)
ANNOVAR_OPERATIONS: str = "g,f,f,f,f,f,f,f"

# Hard-filter expressions for SNPs (GATK best-practice)
SNP_FILTERS: list[tuple[str, str]] = [
    ("QD < 2.0",             "SNP_LowQD"),
    ("FS > 60.0",            "SNP_StrandBias"),
    ("SOR > 3.0",            "SNP_StrandOddsRatio"),
    ("MQ < 40.0",            "SNP_LowMQ"),
    ("MQRankSum < -12.5",    "SNP_MQRankSum"),
    ("ReadPosRankSum < -8.0","SNP_ReadPosRankSum"),
    ("DP < 10",              "LowDepth"),
]
SNP_GENOTYPE_FILTERS: list[tuple[str, str]] = [
    ("GQ < 20", "LowGQ"),
]

# Hard-filter expressions for INDELs (GATK best-practice)
INDEL_FILTERS: list[tuple[str, str]] = [
    ("QD < 2.0",               "INDEL_LowQD"),
    ("FS > 200.0",             "INDEL_StrandBias"),
    ("SOR > 10.0",             "INDEL_StrandOddsRatio"),
    ("ReadPosRankSum < -20.0", "INDEL_ReadPosRankSum"),
    ("DP < 10",                "LowDepth"),
]
INDEL_GENOTYPE_FILTERS: list[tuple[str, str]] = [
    ("GQ < 20", "LowGQ"),
]


@dataclass
class Config:
    """
    Central configuration object for the ExomeFlow pipeline.

    Every field maps to one (or more) CLI arguments.  Sensible defaults are
    provided where the Bash script had hard-coded values.
    """

    # ------------------------------------------------------------------ paths
    input_dir: Path = Path("fastq")
    output_dir: Path = Path("results")

    # Reference files
    reference: Path = Path("hg38.fa")
    dbsnp: Path = Path("dbsnp.vcf.gz")
    mills: Path = Path("Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
    known_indels: Path = Path("Homo_sapiens_assembly38.known_indels.vcf.gz")
    intervals: Path = Path("")           # empty → whole-genome mode

    # ANNOVAR
    annovar_bin: Path = Path("/usr/local/bin/annovar")
    annovar_db: Path = Path("/annovar/humandb")
    annovar_protocols: str = ANNOVAR_PROTOCOLS
    annovar_operations: str = ANNOVAR_OPERATIONS

    # ----------------------------------------------------------- performance
    threads: int = 24
    fastp_threads: int = 8
    annovar_threads: int = 24
    max_workers: int = 1          # parallel samples
    interval_padding: int = 100
    java_opts: str = "-Xmx80g"

    # ---------------------------------------------------------- derived dirs
    # These are set automatically in __post_init__ from output_dir
    qc_dir: Path = field(init=False)
    fastp_dir: Path = field(init=False)
    map_dir: Path = field(init=False)
    vcf_dir: Path = field(init=False)
    log_dir: Path = field(init=False)
    checkpoint_dir: Path = field(init=False)

    def __post_init__(self) -> None:
        self.output_dir = Path(self.output_dir)
        self.qc_dir = self.output_dir / "QC"
        self.fastp_dir = self.output_dir / "filtered_fastp"
        self.map_dir = self.output_dir / "Mapsam"
        self.vcf_dir = self.output_dir / "VCF"
        self.log_dir = self.output_dir / "logs"
        self.checkpoint_dir = self.output_dir / ".checkpoints"

    # -------------------------------------------------------------- helpers
    def setup_directories(self) -> None:
        """Create all output directories that do not yet exist."""
        for d in (
            self.output_dir,
            self.qc_dir,
            self.fastp_dir,
            self.map_dir,
            self.vcf_dir,
            self.log_dir,
            self.checkpoint_dir,
        ):
            d.mkdir(parents=True, exist_ok=True)

    def env(self) -> dict[str, str]:
        """Return an environment dict suitable for subprocess calls."""
        env = os.environ.copy()
        env["JAVA_OPTS"] = self.java_opts
        return env
