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
