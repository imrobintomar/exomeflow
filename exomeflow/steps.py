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
    # Return value matters: an explicit `False` (as opposed to the implicit
    # `None` from a step whose contract is "raise on failure, otherwise
    # succeeded") tells the caller a graceful skip happened and the step
    # must not be checkpointed as done.
    run: Callable[[list[str], "Config"], bool | None]
    applies: Callable[["Config"], bool] = field(default=_always)
