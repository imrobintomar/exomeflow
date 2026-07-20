"""
CLI argument-validation tests. These only exercise Typer/Click's own option
parsing (exists=True path checks) — nothing here reaches run_command's body,
so no real dependency check, download, or pipeline execution ever happens.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

from exomeflow.cli import app

runner = CliRunner()


@pytest.fixture
def valid_input_dir(tmp_path: Path) -> Path:
    d = tmp_path / "fastq"
    d.mkdir()
    (d / "s1_1.fastq.gz").touch()
    (d / "s1_2.fastq.gz").touch()
    return d


@pytest.mark.parametrize(
    "flag", ["--reference", "--dbsnp", "--mills", "--known-indels",
             "--germline-resource", "--panel-of-normals"]
)
def test_run_rejects_nonexistent_file_override(valid_input_dir: Path, tmp_path: Path, flag: str):
    """
    Regression test: found via audit — these overrides used to have no
    existence check at all, so a typo'd path bypassed the dependency table
    (which only ever inspects saved config, never a CLI override) and only
    failed deep inside BWA MEM/GATK instead of failing fast like --input-dir.
    """
    result = runner.invoke(app, [
        "run", "--input-dir", str(valid_input_dir),
        flag, str(tmp_path / "does_not_exist.vcf"),
    ])
    assert result.exit_code == 2
    assert "does not exist" in result.output


@pytest.mark.parametrize("flag", ["--annovar-bin", "--annovar-db"])
def test_run_rejects_nonexistent_directory_override(valid_input_dir: Path, tmp_path: Path, flag: str):
    result = runner.invoke(app, [
        "run", "--input-dir", str(valid_input_dir),
        flag, str(tmp_path / "does_not_exist_dir"),
    ])
    assert result.exit_code == 2
    assert "does not exist" in result.output


def test_run_rejects_nonexistent_intervals(valid_input_dir: Path, tmp_path: Path):
    """
    Regression test: found via audit — a typo'd --intervals path used to
    pass through unvalidated and silently degrade the whole run from WES to
    whole-genome mode (Config.has_intervals treats "given but missing" the
    same as "not given"), with no error or warning anywhere.
    """
    result = runner.invoke(app, [
        "run", "--input-dir", str(valid_input_dir),
        "--intervals", str(tmp_path / "typo.bed"),
    ])
    assert result.exit_code == 2
    assert "does not exist" in result.output


def test_run_accepts_existing_file_override_past_argument_parsing(valid_input_dir: Path, tmp_path: Path):
    """A real, existing path must not be rejected at the argument-parsing
    stage — confirms exists=True doesn't also block legitimate values."""
    real_vcf = tmp_path / "real.vcf.gz"
    real_vcf.touch()
    result = runner.invoke(app, [
        "run", "--input-dir", str(valid_input_dir), "--reference", str(real_vcf),
    ])
    # It should get past Click's own validation — any failure from here on
    # is run_command's own logic (dependency checks etc.), not a rejected
    # CLI argument, so exit code 2 (Click's usage-error code) must not occur.
    assert result.exit_code != 2 or "does not exist" not in result.output


def test_run_rejects_output_path_that_is_an_existing_file(valid_input_dir: Path, tmp_path: Path):
    """
    Regression test: found via audit — unlike --input-dir, --output can't
    use Typer's exists=True (it's documented as "created if absent"), but a
    value that DOES exist and isn't a directory used to pass straight
    through and only fail deep inside Config.setup_directories()'s mkdir()
    call — after the potentially hours-long first-run setup wizard had
    already completed.
    """
    output_path = tmp_path / "not_a_directory"
    output_path.write_text("oops, this is a file")

    result = runner.invoke(app, [
        "run", "--input-dir", str(valid_input_dir), "--output", str(output_path),
    ])
    assert result.exit_code == 1
    assert "not a directory" in result.output
