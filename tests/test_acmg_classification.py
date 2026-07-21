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


# Real InterVar output, confirmed live against an actual run: the header's
# first column is "#Chr" (leading '#', unlike the main table's bare "Chr"),
# and the classification column's raw value is
# " InterVar: <classification words> PVS1=... PS=[...] ..." — the label
# "InterVar:" is IN the cell value, not just the header name.
_REAL_INTERVAR_HEADER = (
    "#Chr\tStart\tEnd\tRef\tAlt\t InterVar: InterVar and Evidence \n"
)


def test_merge_acmg_skips_when_enriched_file_missing(tmp_path: Path, capfd):
    intervar_table = tmp_path / "out.hg38_multianno.txt.intervar"
    intervar_table.write_text(
        _REAL_INTERVAR_HEADER + "1\t100\t100\tA\tG\t InterVar: Pathogenic PVS1=1 PS1=0 \n"
    )
    assert mod._merge_acmg("s1", intervar_table, tmp_path / "does_not_exist.hpo.txt") is False


def test_merge_acmg_appends_classification_columns(tmp_path: Path):
    """
    Regression test: found live — a first-whitespace split used to grab the
    literal label "InterVar:" as the classification, leaving the real
    (often multi-word, e.g. "Likely benign") text stuck inside "evidence"
    instead. Every row's ACMG_classification used to be the same useless
    string "InterVar:".
    """
    intervar_table = tmp_path / "out.hg38_multianno.txt.intervar"
    intervar_table.write_text(
        _REAL_INTERVAR_HEADER + "1\t100\t100\tA\tG\t InterVar: Likely pathogenic PVS1=1 PS1=0 \n"
    )
    enriched = tmp_path / "s1.annovar.hpo.txt"
    enriched.write_text("Chr\tStart\tEnd\tRef\tAlt\n1\t100\t100\tA\tG\n")

    assert mod._merge_acmg("s1", intervar_table, enriched) is True

    lines = enriched.read_text().splitlines()
    header = lines[0].split("\t")
    assert "ACMG_classification" in header and "ACMG_evidence" in header
    row = lines[1].split("\t")
    assert row[header.index("ACMG_classification")] == "Likely pathogenic"
    assert row[header.index("ACMG_evidence")] == "PVS1=1 PS1=0"


def test_merge_acmg_returns_false_when_classification_column_not_found(tmp_path: Path):
    intervar_table = tmp_path / "out.hg38_multianno.txt.intervar"
    intervar_table.write_text("#Chr\tStart\tEnd\tRef\tAlt\tSomeOtherColumn\n1\t100\t100\tA\tG\tx\n")
    enriched = tmp_path / "s1.annovar.hpo.txt"
    enriched.write_text("Chr\tStart\tEnd\tRef\tAlt\n1\t100\t100\tA\tG\n")

    assert mod._merge_acmg("s1", intervar_table, enriched) is False


def test_merge_acmg_uses_hash_chr_and_deduplicates_transcript_rows(tmp_path: Path):
    """
    Regression test: found live against a real 937-variant somatic run —
    InterVar's own output has one row per overlapping transcript/gene per
    variant (standard ANNOVAR gene-based annotation behavior), its first
    column is "#Chr" not "Chr", and its Chr values are bare numbers ("1",
    "2", ...) while the main annotated table uses UCSC-style "chr"-prefixed
    values ("chr1", "chr2", ...). Getting any one of these three wrong
    either drops chromosome from the join key entirely, or makes the join
    match almost nothing (observed live: only 4 of 937 real variants
    matched, all coincidentally-unprefixed HLA allele contigs) — and
    merging without deduplicating first multiplies rows 1:N instead of 1:1
    (937 real variants became 1261 output rows in the same live run).
    """
    intervar_table = tmp_path / "out.hg38_multianno.txt.intervar"
    intervar_table.write_text(
        _REAL_INTERVAR_HEADER
        + "1\t100\t100\tA\tG\t InterVar: Benign PVS1=0 PS1=0 \n"
        + "1\t100\t100\tA\tG\t InterVar: Benign PVS1=0 PS1=0 \n"  # 2nd transcript, same variant
        + "2\t100\t100\tA\tG\t InterVar: Pathogenic PVS1=1 PS1=1 \n"  # different chromosome, same pos
    )
    enriched = tmp_path / "s1.annovar.hpo.txt"
    enriched.write_text(
        "Chr\tStart\tEnd\tRef\tAlt\n"
        "chr1\t100\t100\tA\tG\n"
        "chr2\t100\t100\tA\tG\n"
    )

    assert mod._merge_acmg("s1", intervar_table, enriched) is True

    lines = enriched.read_text().splitlines()
    assert len(lines) == 3  # header + exactly 2 variant rows, not 3
    header = lines[0].split("\t")
    assert "_chr_norm" not in header  # internal join key must not leak into output
    chr1_row = next(r for r in lines[1:] if r.split("\t")[header.index("Chr")] == "chr1")
    chr2_row = next(r for r in lines[1:] if r.split("\t")[header.index("Chr")] == "chr2")
    assert chr1_row.split("\t")[header.index("ACMG_classification")] == "Benign"
    assert chr2_row.split("\t")[header.index("ACMG_classification")] == "Pathogenic"


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


def test_resolve_intervar_db_prefers_own_humandb_when_more_complete(tmp_path: Path, cfg):
    """
    Regression test: found live — a machine with a legacy standalone
    InterVar install can have intervar_bin/humandb already fully populated
    with InterVar's own required (older) database versions, while the
    shared --annovar-db only has the main pipeline's newer versions.
    Unconditionally preferring the shared dir meant InterVar ignored an
    already-complete local set and tried to re-download a ~48GB file from
    ANNOVAR's own (observed to be extremely slow) server instead.
    """
    intervar_bin = tmp_path / "intervar"
    own_db = intervar_bin / "humandb"
    own_db.mkdir(parents=True)
    (intervar_bin / "config.ini").write_text(
        "[Annovar]\ndatabase_names = refGene dbnsfp42a avsnp147\n"
    )
    (own_db / "hg38_refGene.txt").touch()
    (own_db / "hg38_dbnsfp42a.txt").touch()
    (own_db / "hg38_avsnp147.txt").touch()

    shared_db = tmp_path / "shared_annovar_db"
    shared_db.mkdir()
    (shared_db / "hg38_refGene.txt").touch()  # only 1 of 3 required files

    cfg.annovar_db = shared_db
    assert mod._resolve_intervar_db(intervar_bin, cfg) == own_db


def test_resolve_intervar_db_falls_back_to_shared_when_own_is_empty(tmp_path: Path, cfg):
    intervar_bin = tmp_path / "intervar"
    (intervar_bin / "humandb").mkdir(parents=True)
    (intervar_bin / "config.ini").write_text("[Annovar]\ndatabase_names = refGene\n")

    shared_db = tmp_path / "shared_annovar_db"
    shared_db.mkdir()
    cfg.annovar_db = shared_db

    assert mod._resolve_intervar_db(intervar_bin, cfg) == shared_db
