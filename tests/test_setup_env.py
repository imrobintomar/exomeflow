import json
from pathlib import Path

import exomeflow.setup_env as setup_env
from exomeflow.setup_env import (
    ANNOVAR_DATABASES_BY_BUILD,
    annovar_databases_complete,
    detect_annovar_bin,
    detect_gatk_bin,
)


def test_annovar_databases_complete_detects_missing_file(tmp_path: Path):
    db_dir = tmp_path / "humandb"
    db_dir.mkdir()
    for name, _, _ in ANNOVAR_DATABASES_BY_BUILD["hg38"]:
        (db_dir / f"hg38_{name}.txt").touch()

    complete, missing = annovar_databases_complete(db_dir, "hg38")
    assert complete
    assert missing == []

    (db_dir / "hg38_clinvar_20240611.txt").unlink()
    complete, missing = annovar_databases_complete(db_dir, "hg38")
    assert not complete
    assert "clinvar_20240611" in missing


def test_annovar_databases_complete_is_build_aware(tmp_path: Path):
    db_dir = tmp_path / "humandb"
    db_dir.mkdir()
    for name, _, _ in ANNOVAR_DATABASES_BY_BUILD["hg38"]:
        (db_dir / f"hg38_{name}.txt").touch()

    # hg38 files present, but requesting GRCh37 completeness must not be
    # satisfied by them (different buildver prefix, different db names).
    complete, missing = annovar_databases_complete(db_dir, "GRCh37")
    assert not complete
    assert "gnomad211_exome" in missing


def _write_config(config_path: Path, data: dict) -> None:
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(json.dumps(data))


def test_detect_annovar_bin_uses_saved_config_regardless_of_cwd(tmp_path, monkeypatch):
    annovar_dir = tmp_path / "somewhere" / "annovar"
    annovar_dir.mkdir(parents=True)
    (annovar_dir / "table_annovar.pl").touch()

    config_path = tmp_path / "config.json"
    _write_config(config_path, {"annovar_bin": str(annovar_dir)})
    monkeypatch.setattr(setup_env, "CONFIG_PATH", config_path)

    # A directory with no local "annovar" folder — the cwd/home/opt
    # heuristics all miss, so only the saved-config check can succeed.
    other_dir = tmp_path / "unrelated_cwd"
    other_dir.mkdir()
    monkeypatch.chdir(other_dir)
    monkeypatch.setattr(Path, "home", lambda: other_dir)

    assert detect_annovar_bin() == annovar_dir


def test_detect_annovar_bin_ignores_stale_config_path(tmp_path, monkeypatch):
    stale_path = tmp_path / "deleted_annovar"
    config_path = tmp_path / "config.json"
    _write_config(config_path, {"annovar_bin": str(stale_path)})
    monkeypatch.setattr(setup_env, "CONFIG_PATH", config_path)

    other_dir = tmp_path / "unrelated_cwd"
    other_dir.mkdir()
    monkeypatch.chdir(other_dir)
    monkeypatch.setattr(Path, "home", lambda: other_dir)

    # A saved path that no longer exists on disk must never be trusted —
    # whatever comes back (including None) must not be the stale path.
    assert detect_annovar_bin() != stale_path


def test_detect_gatk_bin_uses_saved_config_regardless_of_cwd(tmp_path, monkeypatch):
    gatk_bin = tmp_path / "somewhere" / "gatk-4.6.2.0" / "gatk"
    gatk_bin.parent.mkdir(parents=True)
    gatk_bin.touch()

    config_path = tmp_path / "config.json"
    _write_config(config_path, {"gatk_bin": str(gatk_bin)})
    monkeypatch.setattr(setup_env, "CONFIG_PATH", config_path)

    other_dir = tmp_path / "unrelated_cwd"
    other_dir.mkdir()
    monkeypatch.chdir(other_dir)
    monkeypatch.setattr(Path, "home", lambda: other_dir)
    monkeypatch.setattr(setup_env.shutil, "which", lambda name: None)

    assert detect_gatk_bin() == gatk_bin
