import json
from pathlib import Path

import exomeflow.setup_env as setup_env
from exomeflow.setup_env import (
    ANNOVAR_DATABASES_BY_BUILD,
    _ensure_mim2gene,
    _step_annovar_databases,
    _step_bootstrap_micromamba,
    _step_bundled_tools,
    _step_system_tools,
    annovar_databases_complete,
    detect_annovar_bin,
    detect_annovar_humandb,
    detect_gatk_bin,
    run_setup,
)


def _touch_complete_humandb(db_dir: Path, genome_build: str) -> None:
    """Create .txt (+ .idx for everything but refGene) for every required db."""
    buildver = "hg19" if genome_build == "GRCh37" else "hg38"
    for name, _, _ in ANNOVAR_DATABASES_BY_BUILD[genome_build]:
        txt = db_dir / f"{buildver}_{name}.txt"
        txt.touch()
        if name != "refGene":
            Path(str(txt) + ".idx").touch()


def test_annovar_databases_complete_detects_missing_file(tmp_path: Path):
    db_dir = tmp_path / "humandb"
    db_dir.mkdir()
    _touch_complete_humandb(db_dir, "hg38")

    complete, missing = annovar_databases_complete(db_dir, "hg38")
    assert complete
    assert missing == []

    (db_dir / "hg38_clinvar_20240611.txt").unlink()
    complete, missing = annovar_databases_complete(db_dir, "hg38")
    assert not complete
    assert "clinvar_20240611" in missing


def test_annovar_databases_complete_detects_missing_idx(tmp_path: Path):
    db_dir = tmp_path / "humandb"
    db_dir.mkdir()
    _touch_complete_humandb(db_dir, "hg38")

    # .txt present but its paired .idx missing (e.g. an interrupted
    # download) must still count as incomplete — refGene is exempt since
    # it's gene-based and has no .idx pair.
    (db_dir / "hg38_avsnp150.txt.idx").unlink()
    complete, missing = annovar_databases_complete(db_dir, "hg38")
    assert not complete
    assert "avsnp150" in missing

    (db_dir / "hg38_refGene.txt.idx").unlink(missing_ok=True)
    complete, missing = annovar_databases_complete(db_dir, "hg38")
    assert "refGene" not in missing


def test_annovar_databases_complete_is_build_aware(tmp_path: Path):
    db_dir = tmp_path / "humandb"
    db_dir.mkdir()
    _touch_complete_humandb(db_dir, "hg38")

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


def test_step_bundled_tools_no_prompt_under_assume_yes(tmp_path, monkeypatch):
    monkeypatch.setattr(setup_env, "CONFIG_PATH", tmp_path / "config.json")
    monkeypatch.setattr(setup_env, "detect_gatk_bin", lambda: None)
    monkeypatch.setattr(setup_env, "_step_gatk_download", lambda: None)
    monkeypatch.setattr(setup_env, "detect_annovar_bin", lambda: None)

    def _fail_if_called(*a, **k):
        raise AssertionError("input() should never be called under assume_yes=True")

    monkeypatch.setattr("builtins.input", _fail_if_called)

    gatk, annovar = _step_bundled_tools(assume_yes=True)
    assert annovar is None


def test_step_bundled_tools_prompts_for_annovar_path(tmp_path, monkeypatch):
    annovar_dir = tmp_path / "somewhere" / "annovar"
    annovar_dir.mkdir(parents=True)
    (annovar_dir / "table_annovar.pl").touch()

    monkeypatch.setattr(setup_env, "CONFIG_PATH", tmp_path / "config.json")
    monkeypatch.setattr(setup_env, "detect_gatk_bin", lambda: None)
    monkeypatch.setattr(setup_env, "_step_gatk_download", lambda: None)
    monkeypatch.setattr(setup_env, "detect_annovar_bin", lambda: None)
    monkeypatch.setattr(setup_env.sys.stdin, "isatty", lambda: True)

    answers = iter(["y", str(annovar_dir)])
    monkeypatch.setattr("builtins.input", lambda *_: next(answers))

    gatk, annovar = _step_bundled_tools(assume_yes=False)
    assert annovar == annovar_dir


def test_bootstrap_micromamba_returns_cached_binary(tmp_path, monkeypatch):
    bin_dir = tmp_path / "conda" / "bin"
    bin_dir.mkdir(parents=True)
    fake_bin = bin_dir / "micromamba"
    fake_bin.touch()
    monkeypatch.setattr(setup_env, "MICROMAMBA_BIN", fake_bin)

    def _fail_if_called(*a, **k):
        raise AssertionError("should not re-download when already bootstrapped")

    monkeypatch.setattr(setup_env, "_download_file", _fail_if_called)

    assert _step_bootstrap_micromamba() == fake_bin


def test_bootstrap_micromamba_returns_none_on_download_failure(tmp_path, monkeypatch):
    monkeypatch.setattr(setup_env, "MICROMAMBA_DIR", tmp_path / "conda")
    monkeypatch.setattr(setup_env, "MICROMAMBA_BIN", tmp_path / "conda" / "bin" / "micromamba")
    monkeypatch.setattr(setup_env, "_download_file", lambda url, dest: False)

    assert _step_bootstrap_micromamba() is None


def test_step_system_tools_bootstraps_micromamba_when_conda_missing(tmp_path, monkeypatch):
    fake_mm = tmp_path / "micromamba"
    fake_mm.touch()
    env_dir = tmp_path / "env"

    monkeypatch.setattr(setup_env.shutil, "which", lambda name: None)
    monkeypatch.setattr(setup_env, "_step_bootstrap_micromamba", lambda: fake_mm)
    monkeypatch.setattr(setup_env, "MICROMAMBA_ENV", env_dir)

    calls = []

    def _fake_run_silent(cmd, timeout_s=900):
        calls.append(cmd)
        return True

    monkeypatch.setattr(setup_env, "_run_silent", _fake_run_silent)

    failures = _step_system_tools()
    assert failures == []
    assert len(calls) == len(setup_env._SYSTEM_TOOLS)
    for cmd in calls:
        assert cmd[0] == str(fake_mm)
        assert cmd[1] == "create"
        assert "-p" in cmd and str(env_dir) in cmd


def test_step_system_tools_reports_failure_when_bootstrap_fails(monkeypatch):
    monkeypatch.setattr(setup_env.shutil, "which", lambda name: None)
    monkeypatch.setattr(setup_env, "_step_bootstrap_micromamba", lambda: None)

    failures = _step_system_tools()
    assert len(failures) == 1
    assert "conda not found" in failures[0]


def test_step_annovar_databases_downloads_missing_db_on_partial_match(tmp_path, monkeypatch):
    """
    Regression test for a live bug: a humandb directory with 6/7 required
    databases was reported as "found" and the pipeline proceeded without
    ever downloading the 7th (clinvar_20240611) — the completeness check
    only required *some* database to be present, not all of them.
    """
    annovar_bin = tmp_path / "annovar"
    annovar_bin.mkdir()
    (annovar_bin / "annotate_variation.pl").touch()

    db_dir = tmp_path / "humandb"
    db_dir.mkdir()
    buildver = "hg38"
    missing_db = "clinvar_20240611"
    _touch_complete_humandb(db_dir, "hg38")
    (db_dir / f"{buildver}_{missing_db}.txt").unlink()
    (db_dir / f"{buildver}_{missing_db}.txt.idx").unlink()

    downloaded = []

    def _fake_run_visible(cmd, timeout_s=21600):
        # cmd: ["perl", annotate_pl, "-buildver", buildver, "-downdb",
        #       "-webfrom", "annovar", db_name, db_dir]
        db_name = cmd[7]
        downloaded.append(db_name)
        (db_dir / f"{buildver}_{db_name}.txt").touch()
        (db_dir / f"{buildver}_{db_name}.txt.idx").touch()
        return True

    monkeypatch.setattr(setup_env, "_run_visible", _fake_run_visible)

    resolved_db, failures = _step_annovar_databases(
        annovar_bin, db_dir, genome_build="hg38", assume_yes=True,
    )

    assert downloaded == [missing_db]
    assert failures == []
    assert resolved_db == db_dir
    complete, still_missing = annovar_databases_complete(db_dir, "hg38")
    assert complete
    assert still_missing == []


def test_step_annovar_databases_no_download_when_already_complete(tmp_path, monkeypatch):
    annovar_bin = tmp_path / "annovar"
    annovar_bin.mkdir()
    (annovar_bin / "annotate_variation.pl").touch()

    db_dir = tmp_path / "humandb"
    db_dir.mkdir()
    _touch_complete_humandb(db_dir, "hg38")

    def _fail_if_called(*a, **k):
        raise AssertionError("should not attempt any download when already complete")

    monkeypatch.setattr(setup_env, "_run_visible", _fail_if_called)

    resolved_db, failures = _step_annovar_databases(
        annovar_bin, db_dir, genome_build="hg38", assume_yes=True,
    )
    assert resolved_db == db_dir
    assert failures == []


def test_step_annovar_databases_no_duplicate_check_when_paths_match(tmp_path, monkeypatch):
    """
    Regression test for a live bug: default_db and annovar_bin/humandb are
    frequently the exact same path (whenever no explicit --annovar-db was
    given), which printed the identical "Found N/7 ... missing: ..." line
    twice in a row.
    """
    annovar_bin = tmp_path / "annovar"
    annovar_bin.mkdir()
    (annovar_bin / "annotate_variation.pl").touch()
    db_dir = annovar_bin / "humandb"
    db_dir.mkdir()

    reports = []
    monkeypatch.setattr(setup_env.console, "print", lambda *a, **k: reports.append(a))
    monkeypatch.setattr(setup_env, "detect_annovar_humandb", lambda buildver: None)
    monkeypatch.setattr(setup_env, "_ask", lambda *a, **k: False)

    _step_annovar_databases(annovar_bin, db_dir, genome_build="hg38", assume_yes=True)

    zero_of_seven = [r for r in reports if any("0/7" in str(x) for x in r)]
    assert len(zero_of_seven) == 1


def test_run_setup_reuses_saved_annovar_db(tmp_path, monkeypatch):
    """
    Regression test for a live bug: a previously-saved annovar_db (answered
    interactively in an earlier run) was ignored by run_setup(), which
    always fell back to <annovar_bin>/humandb and re-triggered the whole
    find-or-ask flow from scratch on every subsequent `exomeflow setup`.
    """
    annovar_bin = tmp_path / "annovar"
    annovar_bin.mkdir()
    saved_db = tmp_path / "real_humandb_elsewhere"
    saved_db.mkdir()

    monkeypatch.setattr(setup_env, "CONFIG_PATH", tmp_path / "config.json")
    setup_env.save_config({"annovar_db": str(saved_db)})

    monkeypatch.setattr(setup_env, "_step_bundled_tools", lambda assume_yes=False: (None, annovar_bin))
    monkeypatch.setattr(setup_env, "_step_system_tools", lambda *a, **k: [])
    monkeypatch.setattr(setup_env, "_step_reference_files", lambda *a, **k: (None, []))
    monkeypatch.setattr(setup_env, "_step_intervar", lambda *a, **k: None)
    monkeypatch.setattr(setup_env, "_step_hpo_mapping", lambda: True)
    monkeypatch.setattr(setup_env, "_step_multiqc", lambda: True)

    seen_default_db = {}

    def _fake_step_annovar_databases(annovar_bin_arg, default_db, genome_build, assume_yes=False):
        seen_default_db["path"] = default_db
        return default_db, []

    monkeypatch.setattr(setup_env, "_step_annovar_databases", _fake_step_annovar_databases)

    run_setup(refs_dir=tmp_path / "refs", genome_build="hg38", assume_yes=True)

    assert seen_default_db["path"] == saved_db


def test_detect_annovar_humandb_prefers_most_complete_hit(tmp_path, monkeypatch):
    """
    Regression test: multiple refGene.txt hits on the system (e.g. InterVar
    bundles its own small humandb subset) used to just take find's first
    output line, which isn't a meaningful ordering.
    """
    small_hit = tmp_path / "intervar_humandb"
    small_hit.mkdir()
    (small_hit / "hg38_refGene.txt").touch()

    full_hit = tmp_path / "real_humandb"
    full_hit.mkdir()
    for name in ["refGene", "avsnp150", "clinvar_20240611", "gnomad41_exome"]:
        (full_hit / f"hg38_{name}.txt").touch()

    monkeypatch.setattr(setup_env.shutil, "which", lambda name: "/usr/bin/find")

    class _FakeResult:
        stdout = f"{small_hit}/hg38_refGene.txt\n{full_hit}/hg38_refGene.txt\n"

    monkeypatch.setattr(setup_env.subprocess, "run", lambda *a, **k: _FakeResult())
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "no_fast_candidates_here")

    result = detect_annovar_humandb("hg38")
    assert result == full_hit


def test_detect_annovar_humandb_salvages_partial_output_on_timeout(tmp_path, monkeypatch):
    hit_dir = tmp_path / "humandb"
    hit_dir.mkdir()
    (hit_dir / "hg38_refGene.txt").touch()

    monkeypatch.setattr(setup_env.shutil, "which", lambda name: "/usr/bin/find")
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "no_fast_candidates_here")

    def _raise_timeout(*a, **k):
        raise setup_env.subprocess.TimeoutExpired(
            cmd="find", timeout=30, output=f"{hit_dir}/hg38_refGene.txt\n"
        )

    monkeypatch.setattr(setup_env.subprocess, "run", _raise_timeout)

    result = detect_annovar_humandb("hg38")
    assert result == hit_dir


def test_ensure_mim2gene_downloads_when_missing(tmp_path, monkeypatch):
    """
    Regression test for a live bug: InterVar hard-requires
    intervardb/mim2gene.txt to run its ACMG classification at all, but a
    plain `git clone` of InterVar's repo never provisions it (InterVar's
    own repo doesn't ship it), so classification silently failed on every
    sample with no indication why. Unlike ANNOVAR, this file needs no
    registration — OMIM publishes it as a plain public download.
    """
    intervar_dir = tmp_path / "intervar"
    intervar_dir.mkdir()

    downloaded = {}

    def _fake_download(url, dest):
        downloaded["url"] = url
        downloaded["dest"] = dest
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text("dummy")
        return True

    monkeypatch.setattr(setup_env, "_download_file", _fake_download)

    _ensure_mim2gene(intervar_dir)

    assert downloaded["dest"] == intervar_dir / "intervardb" / "mim2gene.txt"
    assert downloaded["url"] == setup_env.MIM2GENE_URL
    assert (intervar_dir / "intervardb" / "mim2gene.txt").exists()


def test_ensure_mim2gene_skips_when_already_present(tmp_path, monkeypatch):
    intervar_dir = tmp_path / "intervar"
    db_dir = intervar_dir / "intervardb"
    db_dir.mkdir(parents=True)
    (db_dir / "mim2gene.txt").write_text("already here")

    def _fail_if_called(*a, **k):
        raise AssertionError("should not re-download when already present")

    monkeypatch.setattr(setup_env, "_download_file", _fail_if_called)

    _ensure_mim2gene(intervar_dir)
