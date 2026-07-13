from pathlib import Path

from exomeflow.setup_env import ANNOVAR_DATABASES_BY_BUILD, annovar_databases_complete


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
