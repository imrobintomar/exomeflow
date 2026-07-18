from pathlib import Path

from exomeflow.utils import SystemResources, detect_system_resources, recommend_java_opts, recommend_threads


def test_detect_system_resources_handles_nonexistent_output_dir(tmp_path: Path):
    # --output is "created if absent" -- disk_path won't exist yet on a fresh
    # run. Regression test for a real crash found during smoke testing:
    # shutil.disk_usage() raises FileNotFoundError on a missing path.
    missing = tmp_path / "does" / "not" / "exist" / "yet"
    resources = detect_system_resources(missing)
    assert resources.free_disk_gb > 0
    assert resources.cpu_count >= 1


def test_recommend_threads_leaves_headroom_and_caps_at_24():
    assert recommend_threads(SystemResources(48, 64, 48, 1000)) == 24
    assert recommend_threads(SystemResources(8, 16, 12, 1000)) == 6
    assert recommend_threads(SystemResources(2, 8, 6, 1000)) == 1  # never below 1
    assert recommend_threads(SystemResources(1, 4, 2, 1000)) == 1


def test_recommend_java_opts_uses_available_ram_with_floor_and_cap():
    assert recommend_java_opts(SystemResources(48, 64, 47, 1000)) == "-Xmx28g"
    assert recommend_java_opts(SystemResources(4, 8, 5, 1000)) == "-Xmx4g"  # floor
    assert recommend_java_opts(SystemResources(64, 512, 400, 1000)) == "-Xmx80g"  # cap
    assert recommend_java_opts(SystemResources(4, 0, 0, 1000)) == "-Xmx8g"  # unknown RAM


def test_recommend_threads_divides_by_max_workers():
    """
    Regression test: found via audit — thread/heap sizing used to be
    computed for a single worker's exclusive use of the whole machine, with
    no awareness of --max-workers, so N parallel samples each launched
    their own full-machine-sized BWA/GATK processes (Nx oversubscription).
    """
    # (48-2) cores split across 4 workers = 11 each, well under the 24 cap.
    assert recommend_threads(SystemResources(48, 64, 48, 1000), max_workers=4) == 11
    # Per-worker cap still applies even with plenty of cores to go around.
    assert recommend_threads(SystemResources(200, 256, 200, 1000), max_workers=4) == 24
    # Never below 1 per worker even when workers outnumber available cores.
    assert recommend_threads(SystemResources(4, 8, 4, 1000), max_workers=8) == 1
    # max_workers=1 (the default) matches the pre-existing single-worker behavior.
    assert recommend_threads(SystemResources(48, 64, 48, 1000), max_workers=1) == 24


def test_recommend_java_opts_divides_by_max_workers():
    assert recommend_java_opts(SystemResources(48, 64, 48, 1000), max_workers=4) == "-Xmx7g"
    # Floor still applies per worker, not just in aggregate.
    assert recommend_java_opts(SystemResources(48, 64, 8, 1000), max_workers=4) == "-Xmx4g"
    assert recommend_java_opts(SystemResources(48, 64, 48, 1000), max_workers=1) == "-Xmx28g"
