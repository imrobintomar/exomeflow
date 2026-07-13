"""
Pluggable pipeline step registry.

A SampleStep runs once per sample; a CohortStep runs once after every
sample's SampleSteps have finished. pipeline.py builds an explicit ordered
list of each from the step functions in the other modules, gating each one
with an `applies(cfg)` predicate instead of hardcoding branches per mode.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Callable

if TYPE_CHECKING:
    from exomeflow.config import Config
    from exomeflow.utils import Checkpoint


def _always(cfg: "Config") -> bool:
    return True


@dataclass(frozen=True)
class SampleStep:
    name: str
    run: Callable[[str, "Config", "Checkpoint"], None]
    applies: Callable[["Config"], bool] = field(default=_always)


@dataclass(frozen=True)
class CohortStep:
    name: str
    run: Callable[[list[str], "Config"], None]
    applies: Callable[["Config"], bool] = field(default=_always)
