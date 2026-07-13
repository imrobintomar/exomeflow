from pathlib import Path

import exomeflow.acmg_classification as mod


def test_run_intervar_tool_skips_gracefully_when_not_installed(tmp_path: Path, cfg, monkeypatch):
    monkeypatch.setattr("exomeflow.setup_env.detect_intervar_bin", lambda: None)
    result = mod._run_intervar_tool("s1", tmp_path / "s1_PASS.vcf", tmp_path / "s1.intervar", cfg)
    assert result is None


def test_merge_acmg_skips_when_enriched_file_missing(tmp_path: Path, capfd):
    intervar_table = tmp_path / "out.hg38_multianno.txt.intervar"
    intervar_table.write_text(
        "Chr\tStart\tEnd\tRef\tAlt\t InterVar: InterVar and Evidence \n"
        "1\t100\t100\tA\tG\tPathogenic PVS1=1 PS1=0\n"
    )
    mod._merge_acmg("s1", intervar_table, tmp_path / "does_not_exist.hpo.txt")
    # must not raise; nothing to assert beyond "no exception"


def test_merge_acmg_appends_classification_columns(tmp_path: Path):
    intervar_table = tmp_path / "out.hg38_multianno.txt.intervar"
    intervar_table.write_text(
        "Chr\tStart\tEnd\tRef\tAlt\t InterVar: InterVar and Evidence \n"
        "1\t100\t100\tA\tG\tPathogenic PVS1=1 PS1=0\n"
    )
    enriched = tmp_path / "s1.annovar.hpo.txt"
    enriched.write_text("Chr\tStart\tEnd\tRef\tAlt\n1\t100\t100\tA\tG\n")

    mod._merge_acmg("s1", intervar_table, enriched)

    lines = enriched.read_text().splitlines()
    header = lines[0].split("\t")
    assert "ACMG_classification" in header and "ACMG_evidence" in header
    row = lines[1].split("\t")
    assert row[header.index("ACMG_classification")] == "Pathogenic"
    assert row[header.index("ACMG_evidence")] == "PVS1=1 PS1=0"
