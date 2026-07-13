from pathlib import Path

import pytest

from exomeflow.config import Config


@pytest.fixture
def cfg(tmp_path: Path) -> Config:
    return Config(
        input_dir=tmp_path / "fastq",
        output_dir=tmp_path / "results",
        reference=tmp_path / "ref.fa",
        dbsnp=tmp_path / "dbsnp.vcf.gz",
        mills=tmp_path / "mills.vcf.gz",
        known_indels=tmp_path / "known_indels.vcf.gz",
        annovar_bin=tmp_path / "annovar",
        annovar_db=tmp_path / "annovar" / "humandb",
    )
