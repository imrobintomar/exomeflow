import subprocess
from pathlib import Path

import exomeflow.acmg_classification as mod
from exomeflow.utils import Checkpoint


def test_run_intervar_tool_skips_gracefully_when_not_installed(tmp_path: Path, cfg, monkeypatch):
    monkeypatch.setattr("exomeflow.setup_env.detect_intervar_bin", lambda: None)
    result = mod._run_intervar_tool("s1", tmp_path / "s1_PASS.vcf", tmp_path / "s1.intervar", cfg)
    assert result is None


def test_run_intervar_tool_shares_annovar_db_and_passes_intervardb_flag(tmp_path: Path, cfg, monkeypatch):
    """
    Regression test for two live bugs: (1) InterVar was pointed at its own
    fully separate humandb copy instead of reusing cfg.annovar_db, wasting
    tens of GB on files (refGene, ensGene, knownGene) already downloaded by
    the main pipeline; (2) -t/--database_intervar (InterVar's own ACMG
    criteria database, distinct from -d) was never passed at all, so
    InterVar fell back to a relative default path resolved against the
    pipeline's cwd instead of its own install directory.
    """
    intervar_bin = tmp_path / "intervar"
    (intervar_bin / "intervardb").mkdir(parents=True)
    (intervar_bin / "Intervar.py").touch()
    monkeypatch.setattr("exomeflow.setup_env.detect_intervar_bin", lambda: intervar_bin)

    captured = {}

    class _FakeProc:
        returncode = 0

        def communicate(self, timeout=None):
            return "", ""

    def _fake_popen(cmd, **kwargs):
        captured["cmd"] = cmd
        return _FakeProc()

    monkeypatch.setattr(subprocess, "Popen", _fake_popen)

    mod._run_intervar_tool("s1", tmp_path / "s1_PASS.vcf", tmp_path / "s1.intervar", cfg)

    cmd = captured["cmd"]
    assert "-d" in cmd
    assert cmd[cmd.index("-d") + 1] == str(cfg.annovar_db)
    assert "-t" in cmd
    assert cmd[cmd.index("-t") + 1] == str(intervar_bin / "intervardb")


def test_merge_acmg_skips_when_enriched_file_missing(tmp_path: Path, capfd):
    intervar_table = tmp_path / "out.hg38_multianno.txt.intervar"
    intervar_table.write_text(
        "Chr\tStart\tEnd\tRef\tAlt\t InterVar: InterVar and Evidence \n"
        "1\t100\t100\tA\tG\tPathogenic PVS1=1 PS1=0\n"
    )
    assert mod._merge_acmg("s1", intervar_table, tmp_path / "does_not_exist.hpo.txt") is False


def test_merge_acmg_appends_classification_columns(tmp_path: Path):
    intervar_table = tmp_path / "out.hg38_multianno.txt.intervar"
    intervar_table.write_text(
        "Chr\tStart\tEnd\tRef\tAlt\t InterVar: InterVar and Evidence \n"
        "1\t100\t100\tA\tG\tPathogenic PVS1=1 PS1=0\n"
    )
    enriched = tmp_path / "s1.annovar.hpo.txt"
    enriched.write_text("Chr\tStart\tEnd\tRef\tAlt\n1\t100\t100\tA\tG\n")

    assert mod._merge_acmg("s1", intervar_table, enriched) is True

    lines = enriched.read_text().splitlines()
    header = lines[0].split("\t")
    assert "ACMG_classification" in header and "ACMG_evidence" in header
    row = lines[1].split("\t")
    assert row[header.index("ACMG_classification")] == "Pathogenic"
    assert row[header.index("ACMG_evidence")] == "PVS1=1 PS1=0"


def test_merge_acmg_returns_false_when_classification_column_not_found(tmp_path: Path):
    intervar_table = tmp_path / "out.hg38_multianno.txt.intervar"
    intervar_table.write_text("Chr\tStart\tEnd\tRef\tAlt\tSomeOtherColumn\n1\t100\t100\tA\tG\tx\n")
    enriched = tmp_path / "s1.annovar.hpo.txt"
    enriched.write_text("Chr\tStart\tEnd\tRef\tAlt\n1\t100\t100\tA\tG\n")

    assert mod._merge_acmg("s1", intervar_table, enriched) is False


def test_run_intervar_does_not_checkpoint_when_merge_fails(tmp_path: Path, cfg, monkeypatch):
    """
    Regression test: found via audit — run_intervar() used to checkpoint
    unconditionally whenever InterVar produced *any* output file, even if
    _merge_acmg() then failed to actually attach ACMG columns (e.g. its
    classification column couldn't be located). That silently and
    permanently blocked any future retry once the real cause was fixed.
    """
    cfg.vcf_dir.mkdir(parents=True)
    (cfg.vcf_dir / "s1_PASS.vcf").write_text("##fileformat=VCFv4.2\n")

    fake_table = tmp_path / "s1.intervar"
    fake_table.touch()
    monkeypatch.setattr(mod, "_run_intervar_tool", lambda *a, **k: fake_table)
    monkeypatch.setattr(mod, "_merge_acmg", lambda *a, **k: False)

    checkpoint = Checkpoint(cfg.checkpoint_dir)
    mod.run_intervar("s1", cfg, checkpoint)

    assert not checkpoint.done("s1", mod.STEP)


def test_run_intervar_cohort_returns_false_when_intervar_missing(tmp_path: Path, cfg, monkeypatch):
    """
    Regression test: found via audit — run_intervar_cohort() always
    implicitly returned None, which pipeline.py's cohort loop treats as
    success. So with --joint-genotyping, cohort ACMG got checkpointed done
    even when InterVar was missing/unprovisioned for that run.
    """
    monkeypatch.setattr(mod, "_run_intervar_tool", lambda *a, **k: None)
    assert mod.run_intervar_cohort(["s1"], cfg) is False
