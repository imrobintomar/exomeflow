import gzip
from pathlib import Path

import pytest

from exomeflow.utils import (
    Checkpoint,
    PipelineStepError,
    _parse_version,
    _version_ok,
    count_variants,
    detect_samples,
    resolve_fastq_pair,
    run_cmd,
)


def test_detect_samples_pairs_by_prefix(tmp_path: Path):
    for name in ["s1_1.fastq.gz", "s1_2.fastq.gz", "s2_1.fastq.gz", "s2_2.fastq.gz"]:
        (tmp_path / name).touch()
    assert detect_samples(tmp_path) == ["s1", "s2"]


def test_detect_samples_ignores_non_matching_files(tmp_path: Path):
    (tmp_path / "s1_1.fastq.gz").touch()
    (tmp_path / "s1_2.fastq.gz").touch()
    (tmp_path / "readme.txt").touch()
    assert detect_samples(tmp_path) == ["s1"]


def test_detect_samples_raises_when_empty(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        detect_samples(tmp_path)


def test_detect_samples_accepts_r1_r2_convention(tmp_path: Path):
    for name in ["s1_R1.fastq.gz", "s1_R2.fastq.gz", "s2_R1.fastq.gz", "s2_R2.fastq.gz"]:
        (tmp_path / name).touch()
    assert detect_samples(tmp_path) == ["s1", "s2"]


def test_detect_samples_mixed_conventions(tmp_path: Path):
    for name in ["s1_1.fastq.gz", "s1_2.fastq.gz", "s2_R1.fastq.gz", "s2_R2.fastq.gz"]:
        (tmp_path / name).touch()
    assert detect_samples(tmp_path) == ["s1", "s2"]


def test_resolve_fastq_pair_prefers_underscore_numeric(tmp_path: Path):
    (tmp_path / "s1_1.fastq.gz").touch()
    (tmp_path / "s1_2.fastq.gz").touch()
    r1, r2 = resolve_fastq_pair(tmp_path, "s1")
    assert r1.name == "s1_1.fastq.gz"
    assert r2.name == "s1_2.fastq.gz"


def test_resolve_fastq_pair_falls_back_to_r1_r2(tmp_path: Path):
    (tmp_path / "s1_R1.fastq.gz").touch()
    (tmp_path / "s1_R2.fastq.gz").touch()
    r1, r2 = resolve_fastq_pair(tmp_path, "s1")
    assert r1.name == "s1_R1.fastq.gz"
    assert r2.name == "s1_R2.fastq.gz"


def test_resolve_fastq_pair_raises_when_missing(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        resolve_fastq_pair(tmp_path, "nope")


def test_checkpoint_mark_and_done(tmp_path: Path):
    cp = Checkpoint(tmp_path / ".checkpoints")
    assert not cp.done("s1", "fastp")
    cp.mark("s1", "fastp")
    assert cp.done("s1", "fastp")
    assert not cp.done("s2", "fastp")


def test_checkpoint_cohort_namespace(tmp_path: Path):
    cp = Checkpoint(tmp_path / ".checkpoints")
    assert not cp.done("__cohort__", "multiqc")
    cp.mark("__cohort__", "multiqc")
    assert cp.done("__cohort__", "multiqc")


def test_checkpoint_isolates_by_genome_build(tmp_path: Path):
    hg38 = Checkpoint(tmp_path / ".checkpoints", genome_build="hg38")
    grch37 = Checkpoint(tmp_path / ".checkpoints", genome_build="GRCh37")
    hg38.mark("s1", "haplotype")
    assert hg38.done("s1", "haplotype")
    assert not grch37.done("s1", "haplotype")


def test_checkpoint_isolates_by_joint_genotyping(tmp_path: Path):
    """
    Regression test: found via audit — HaplotypeCaller writes a different
    output file (.g.vcf.gz vs .vcf) depending on joint_genotyping, but the
    checkpoint key used to ignore this dimension entirely. Toggling
    --joint-genotyping on an existing --output dir would see "haplotype" as
    already done and never produce the GVCF the cohort phase needs.
    """
    std = Checkpoint(tmp_path / ".checkpoints", joint_genotyping=False)
    jg = Checkpoint(tmp_path / ".checkpoints", joint_genotyping=True)
    std.mark("s1", "haplotype")
    assert std.done("s1", "haplotype")
    assert not jg.done("s1", "haplotype")


def test_checkpoint_isolates_by_mode(tmp_path: Path):
    """
    Regression test: found via audit — germline and somatic filtering both
    write to the same _PASS.vcf filename, but the checkpoint key used to
    ignore mode entirely. Switching --mode on an existing --output dir
    would leave the "annovar" checkpoint from the old mode's run stuck
    "done", so annotation never reran against the new mode's PASS VCF.
    """
    germline = Checkpoint(tmp_path / ".checkpoints", mode="germline")
    somatic = Checkpoint(tmp_path / ".checkpoints", mode="somatic")
    germline.mark("s1", "annovar")
    assert germline.done("s1", "annovar")
    assert not somatic.done("s1", "annovar")


@pytest.mark.parametrize(
    "text,expected",
    [
        ("1.13", (1, 13)),
        ("bwa 0.7.17-r1188", (0, 7, 17)),
        ("no digits here", (0,)),
    ],
)
def test_parse_version(text, expected):
    assert _parse_version(text) == expected


def test_version_ok():
    assert _version_ok(_parse_version("1.13.2"), "1.13")
    assert not _version_ok(_parse_version("1.12"), "1.13")


_VCF_BODY = "##fileformat=VCFv4.2\n#CHROM\tPOS\n" + "chr1\t100\n" * 3


def test_count_variants_plain_text(tmp_path: Path):
    vcf = tmp_path / "s1.vcf"
    vcf.write_text(_VCF_BODY)
    assert count_variants(vcf) == 3


def test_count_variants_bgzipped(tmp_path: Path):
    """
    Regression test: found via audit — GATK auto-bgzips any -O output whose
    filename ends in .gz (e.g. the cohort joint-genotyped VCF), and a plain
    open() on a gzipped file raises UnicodeDecodeError instead of counting.
    """
    vcf = tmp_path / "cohort.vcf.gz"
    with gzip.open(vcf, "wt", encoding="utf-8") as fh:
        fh.write(_VCF_BODY)
    assert count_variants(vcf) == 3


def test_count_variants_missing_file(tmp_path: Path):
    assert count_variants(tmp_path / "does_not_exist.vcf") == 0


def test_run_cmd_raises_on_nonzero_exit():
    """
    Regression test: found via audit — every existing test that touches
    run_cmd() monkeypatches it away entirely, so its actual exit-code check
    (the thing every step function depends on to detect a failed tool) was
    never exercised. A regression here (e.g. an inverted condition) would
    have passed the whole suite undetected.
    """
    import sys
    with pytest.raises(PipelineStepError):
        run_cmd([sys.executable, "-c", "import sys; sys.exit(1)"], step_name="test")


def test_run_cmd_returns_completed_process_on_success():
    import sys
    result = run_cmd([sys.executable, "-c", "pass"], step_name="test")
    assert result.returncode == 0
