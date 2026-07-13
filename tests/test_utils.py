import gzip
from pathlib import Path

import pytest

from exomeflow.utils import Checkpoint, _parse_version, _version_ok, count_variants, detect_samples


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
