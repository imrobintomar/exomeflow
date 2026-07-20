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
from typing import Literal


# ANNOVAR buildver naming differs from GATK/UCSC's for GRCh37 — single
# source of truth (was previously duplicated as an inline ternary in three
# separate places; found via audit).
ANNOVAR_BUILDVER: dict[str, str] = {"hg38": "hg38", "GRCh37": "hg19"}

# ANNOVAR annotation databases, per genome build. Verified against ANNOVAR's
# own `-webfrom annovar` mirror listing (hg38_avdblist.txt / hg19_avdblist.txt)
# at implementation time:
#   - clinvar_20240611 (not _20240416 — that build was never published under
#     that name, so the auto-downloader would silently fail on it).
#   - cosmic84_coding is dropped entirely: COSMIC requires its own separate
#     registered download and was never obtainable via `-webfrom annovar`
#     for either build — add it back manually to --annovar-protocols if you
#     have a licensed COSMIC db.
#   - gnomAD v4.1 (gnomad41_exome/genome) is hg38-only — it was never
#     released for GRCh37/hg19. The hg19 equivalent is gnomAD v2.1.1
#     (gnomad211_exome/genome); using the hg38 names under --genome-build
#     GRCh37 would fail since those files don't exist on that mirror. Found
#     via audit — the shipped default previously used the hg38 names
#     unconditionally, silently breaking GRCh37 annotation.
ANNOVAR_PROTOCOLS_BY_BUILD: dict[str, str] = {
    "hg38": (
        "refGene,dbnsfp47a,clinvar_20240611,gnomad41_exome,"
        "gnomad41_genome,avsnp150,exac03"
    ),
    "GRCh37": (
        "refGene,dbnsfp47a,clinvar_20240611,gnomad211_exome,"
        "gnomad211_genome,avsnp150,exac03"
    ),
}
ANNOVAR_OPERATIONS_BY_BUILD: dict[str, str] = {
    "hg38": "g,f,f,f,f,f,f",
    "GRCh37": "g,f,f,f,f,f,f",
}

# Backward-compatible aliases — Config's dataclass field defaults need a
# concrete value at class-definition time (can't depend on another field);
# __post_init__ swaps these for the GRCh37 variants when genome_build
# requests it and the caller didn't explicitly override the protocol list.
ANNOVAR_PROTOCOLS: str = ANNOVAR_PROTOCOLS_BY_BUILD["hg38"]
ANNOVAR_OPERATIONS: str = ANNOVAR_OPERATIONS_BY_BUILD["hg38"]


def intervals_present(path: Path | None) -> bool:
    """
    True only when *path* is a real, existing file — not None, not an empty
    sentinel. Shared by Config.has_intervals and cli.py's pre-Config
    --joint-genotyping/--cnv validation (previously two independent copies
    of the same predicate; found via audit).
    """
    return path is not None and path.exists()

# Hard-filter expressions for SNPs — exactly GATK's official germline
# short-variant hard-filter recommendation, no added thresholds:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890471
# (A site-level DP<10 filter was previously applied here too, on top of
# GATK's own recommendation. It isn't part of GATK's hard-filter guidance,
# and unlike the other expressions here it doesn't distinguish "modest but
# real" from "artifact" — QD (quality normalized by depth) already covers
# that. Removed so real variants at borderline-but-genuine depth, e.g. near
# capture-kit target edges, aren't discarded by a threshold GATK doesn't
# actually recommend.)
SNP_FILTERS: list[tuple[str, str]] = [
    ("QD < 2.0",             "SNP_LowQD"),
    ("FS > 60.0",            "SNP_StrandBias"),
    ("SOR > 3.0",            "SNP_StrandOddsRatio"),
    ("MQ < 40.0",            "SNP_LowMQ"),
    ("MQRankSum < -12.5",    "SNP_MQRankSum"),
    ("ReadPosRankSum < -8.0", "SNP_ReadPosRankSum"),
]
# Genotype-level (not site-level): flags a sample's own call as low-confidence
# via the VCF's FT field, but does not exclude the variant from *_PASS.vcf —
# SelectVariants --exclude-filtered only checks the site FILTER column.
SNP_GENOTYPE_FILTERS: list[tuple[str, str]] = [
    ("GQ < 20", "LowGQ"),
]

# Hard-filter expressions for INDELs — same GATK official set, see above.
INDEL_FILTERS: list[tuple[str, str]] = [
    ("QD < 2.0",               "INDEL_LowQD"),
    ("FS > 200.0",             "INDEL_StrandBias"),
    ("SOR > 10.0",             "INDEL_StrandOddsRatio"),
    ("ReadPosRankSum < -20.0", "INDEL_ReadPosRankSum"),
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
    # None -> whole-genome mode. NOTE: Path("") normalizes to Path(".") (the
    # cwd), which is always truthy and always .exists() — a falsy *empty*
    # Path can never signal "no intervals given", so this must stay a real
    # None sentinel, not Path(""), or every "was --intervals passed?" check
    # downstream silently short-circuits to "yes, and it's the cwd."
    intervals: Path | None = None

    # ANNOVAR
    annovar_bin: Path = Path("/usr/local/bin/annovar")
    annovar_db: Path = Path("/annovar/humandb")
    annovar_protocols: str = ANNOVAR_PROTOCOLS
    annovar_operations: str = ANNOVAR_OPERATIONS

    # -------------------------------------------------------------- V2 mode
    mode: Literal["germline", "somatic"] = "germline"
    genome_build: Literal["hg38", "GRCh37"] = "hg38"
    joint_genotyping: bool = False
    call_cnv: bool = False
    germline_resource: Path | None = None  # Mutect2 --germline-resource (gnomAD AF-only)
    panel_of_normals: Path | None = None  # Mutect2 --panel-of-normals (pre-built PoN VCF)
    common_biallelic_sites: Path | None = None  # GetPileupSummaries site list (small subset, not germline_resource)

    # ----------------------------------------------------------- performance
    threads: int = 24
    fastp_threads: int = 8
    annovar_threads: int = 24
    max_workers: int = 1          # parallel samples
    interval_padding: int = 100
    java_opts: str = "-Xmx80g"

    # ---------------------------------------------------------- derived dirs
    # These are set automatically in __post_init__ from output_dir
    fastp_dir: Path = field(init=False)
    map_dir: Path = field(init=False)
    vcf_dir: Path = field(init=False)
    cohort_dir: Path = field(init=False)
    cnv_dir: Path = field(init=False)
    log_dir: Path = field(init=False)
    checkpoint_dir: Path = field(init=False)

    def __post_init__(self) -> None:
        # Defense-in-depth, same rationale as the protocols/operations count
        # check below: cli.py already validates these, but any other caller
        # building a Config directly bypassed it entirely. Literal[...] is
        # not enforced at runtime, so a typo'd mode (e.g. "Somatic") used to
        # construct silently, match neither the germline nor somatic step
        # predicates, and produce zero variant calling while every step
        # still reported success. Found via audit.
        if self.mode not in ("germline", "somatic"):
            raise ValueError(f"mode must be 'germline' or 'somatic', got {self.mode!r}")
        if self.genome_build not in ("hg38", "GRCh37"):
            raise ValueError(f"genome_build must be 'hg38' or 'GRCh37', got {self.genome_build!r}")

        self.output_dir = Path(self.output_dir)
        self.fastp_dir = self.output_dir / "filtered_fastp"
        self.map_dir = self.output_dir / "Mapsam"
        self.vcf_dir = self.output_dir / "VCF"
        self.cohort_dir = self.vcf_dir / "cohort"
        self.cnv_dir = self.output_dir / "CNV"
        self.log_dir = self.output_dir / "logs"
        self.checkpoint_dir = self.output_dir / ".checkpoints"

        # GRCh37 needs a different ANNOVAR protocol list (gnomAD v4.1 is
        # hg38-only) — swap only if the caller left the field at its hg38
        # dataclass default; an explicit --annovar-protocols override is
        # left untouched regardless of build. Bug found via audit: this
        # used to be a single hg38-only default applied unconditionally.
        if self.genome_build == "GRCh37":
            if self.annovar_protocols == ANNOVAR_PROTOCOLS_BY_BUILD["hg38"]:
                self.annovar_protocols = ANNOVAR_PROTOCOLS_BY_BUILD["GRCh37"]
            if self.annovar_operations == ANNOVAR_OPERATIONS_BY_BUILD["hg38"]:
                self.annovar_operations = ANNOVAR_OPERATIONS_BY_BUILD["GRCh37"]

        # Defense-in-depth: this was previously validated only in cli.py, so
        # any other caller building a Config directly (tests, a future
        # second CLI command, programmatic use of run_pipeline) bypassed the
        # check entirely and only discovered a mismatch via table_annovar.pl's
        # own cryptic error. Found via audit.
        n_protocols = self.annovar_protocols.count(",") + 1
        n_operations = self.annovar_operations.count(",") + 1
        if n_protocols != n_operations:
            raise ValueError(
                f"annovar_protocols has {n_protocols} entries but "
                f"annovar_operations has {n_operations} — ANNOVAR pairs each "
                f"protocol with an operation, these must match."
            )

    # -------------------------------------------------------------- helpers
    @property
    def has_intervals(self) -> bool:
        """True only when a real --intervals path was given and it exists."""
        return intervals_present(self.intervals)

    @property
    def annovar_buildver(self) -> str:
        """ANNOVAR's own buildver naming differs from the GATK/UCSC one for GRCh37."""
        return ANNOVAR_BUILDVER[self.genome_build]

    def setup_directories(self) -> None:
        """Create all output directories that do not yet exist."""
        dirs = [
            self.output_dir,
            self.fastp_dir,
            self.map_dir,
            self.vcf_dir,
            self.log_dir,
            self.checkpoint_dir,
        ]
        if self.joint_genotyping:
            dirs.append(self.cohort_dir)
        if self.call_cnv:
            dirs.append(self.cnv_dir)
        for d in dirs:
            d.mkdir(parents=True, exist_ok=True)

    def env(self) -> dict[str, str]:
        """Return an environment dict suitable for subprocess calls."""
        env = os.environ.copy()
        env["JAVA_OPTS"] = self.java_opts
        return env
