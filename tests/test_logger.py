import logging
from pathlib import Path

from exomeflow.logger import close_sample_logger, get_sample_logger


def test_close_sample_logger_releases_the_file_handler(tmp_path: Path):
    """
    Regression test for a live-audit finding: a ProcessPoolExecutor worker
    is reused across many samples whenever there are more samples than
    --max-workers, and get_sample_logger opened a fresh FileHandler per
    *new* sample name with nothing ever closing the previous one — a worker
    processing many samples over its lifetime leaked one open file
    descriptor per sample, risking `ulimit -n` exhaustion on large cohorts.
    """
    logger = get_sample_logger("s1", tmp_path, "20260101_000000")
    assert len(logger.handlers) == 1
    handler = logger.handlers[0]
    assert handler.stream is not None and not handler.stream.closed

    close_sample_logger("s1")

    assert logging.getLogger("exomeflow.s1").handlers == []
    assert handler.stream is None or handler.stream.closed


def test_get_sample_logger_reopens_after_close(tmp_path: Path):
    """A sample name can legitimately recur (e.g. across separate
    exomeflow run invocations reusing the same process in tests) — after
    close_sample_logger, get_sample_logger must open a fresh handler rather
    than silently staying handler-less forever."""
    get_sample_logger("s1", tmp_path, "20260101_000000")
    close_sample_logger("s1")

    logger = get_sample_logger("s1", tmp_path, "20260101_000001")
    assert len(logger.handlers) == 1

    close_sample_logger("s1")  # cleanup


def test_close_sample_logger_is_a_noop_for_unknown_sample():
    close_sample_logger("never_existed")  # must not raise
