"""
Structured logging for ExomeFlow.

Provides a pipeline-wide logger (rich console + rotating file) and a
per-sample child logger that writes to its own log file.  A custom
SUCCESS level sits between INFO and WARNING, mirroring the Bash script.
"""

from __future__ import annotations

import logging
from pathlib import Path

from rich.console import Console
from rich.logging import RichHandler

# ---------------------------------------------------------------------------
# Custom SUCCESS level (between INFO=20 and WARNING=30)
# ---------------------------------------------------------------------------
SUCCESS_LEVEL: int = 25
logging.addLevelName(SUCCESS_LEVEL, "SUCCESS")


def _success(self: logging.Logger, message: str, *args, **kwargs) -> None:  # type: ignore[override]
    if self.isEnabledFor(SUCCESS_LEVEL):
        self._log(SUCCESS_LEVEL, message, args, **kwargs)


logging.Logger.success = _success  # type: ignore[attr-defined]

# ---------------------------------------------------------------------------
# Module-level console (shared across all loggers for pretty output)
# ---------------------------------------------------------------------------
_console = Console(stderr=True)


def get_pipeline_logger(log_dir: Path, timestamp: str) -> logging.Logger:
    """
    Return the root pipeline logger.

    Writes to:
      - stderr  (rich, coloured)
      - <log_dir>/analysis_<timestamp>.log  (plain text)
      - <log_dir>/errors_<timestamp>.log    (ERROR+ only)
    """
    logger = logging.getLogger("exomeflow")
    if logger.handlers:          # already configured
        return logger

    logger.setLevel(logging.DEBUG)

    # --- rich stderr handler -----------------------------------------------
    rich_handler = RichHandler(
        console=_console,
        show_time=True,
        show_level=True,
        show_path=False,
        markup=True,
        rich_tracebacks=True,
    )
    rich_handler.setLevel(logging.DEBUG)
    logger.addHandler(rich_handler)

    log_dir.mkdir(parents=True, exist_ok=True)

    # --- plain-text analysis log -------------------------------------------
    fmt = logging.Formatter(
        "[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    analysis_handler = logging.FileHandler(
        log_dir / f"analysis_{timestamp}.log", encoding="utf-8"
    )
    analysis_handler.setLevel(logging.DEBUG)
    analysis_handler.setFormatter(fmt)
    logger.addHandler(analysis_handler)

    # --- error-only log ----------------------------------------------------
    error_handler = logging.FileHandler(
        log_dir / f"errors_{timestamp}.log", encoding="utf-8"
    )
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(fmt)
    logger.addHandler(error_handler)

    return logger


def get_sample_logger(sample: str, log_dir: Path, timestamp: str) -> logging.Logger:
    """
    Return a child logger for a single sample.

    Inherits handlers from the pipeline logger and additionally writes to
    <log_dir>/<sample>_<timestamp>.log.
    """
    logger = logging.getLogger(f"exomeflow.{sample}")
    if logger.handlers:
        return logger

    log_dir.mkdir(parents=True, exist_ok=True)
    fmt = logging.Formatter(
        "[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    sample_handler = logging.FileHandler(
        log_dir / f"{sample}_{timestamp}.log", encoding="utf-8"
    )
    sample_handler.setLevel(logging.DEBUG)
    sample_handler.setFormatter(fmt)
    logger.addHandler(sample_handler)

    return logger
