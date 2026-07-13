"""
ACMG/AMP variant classification via InterVar (github.com/WGLab/InterVar) —
wrapping an established, peer-reviewed classifier rather than re-implementing
ACMG criteria in-house, since misclassifying pathogenicity has real clinical
stakes. Runs automatically after annotation but degrades to a skip-with-
warning if InterVar isn't installed, rather than failing the pipeline.

NOTE: InterVar's exact CLI surface/output columns were reconstructed from its
published interface at implementation time and were not hands-on verified
against a live install (flagged in the V2 plan as an open integration risk).
Verify against your InterVar version and adjust `_INTERVAR_COLUMN_HINT` /
the argv below if the column name or flags differ.
"""

from __future__ import annotations

import logging
import os
import signal
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from exomeflow.utils import Checkpoint

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")

STEP = "acmg"

# Substring used to locate InterVar's classification column, since its header
# includes free text ("InterVar: InterVar and Evidence ") rather than a fixed name.
_INTERVAR_COLUMN_HINT = "InterVar"


_INTERVAR_TIMEOUT_S = 1800  # InterVar's config.ini pins its own DB set, independent
# of --annovar-db/--annovar-protocols, and auto-downloads any it's missing from
# ANNOVAR's server on a first run — bounded here so a slow/stalled download of
# InterVar's own dependencies can never hang the rest of the pipeline.


def _run_intervar_tool(label: str, vcf: Path, out_prefix: Path, cfg: "Config") -> Path | None:
    from exomeflow.setup_env import detect_intervar_bin

    intervar_bin = detect_intervar_bin()
    if intervar_bin is None:
        logger.warning(
            "[%s] InterVar not found — skipping ACMG classification "
            "(auto-installed by `exomeflow run`'s dependency check).",
            label,
        )
        return None

    # InterVar's own database dir, isolated from cfg.annovar_db: InterVar's
    # config.ini requires a fixed set of (specifically versioned) databases of
    # its own, unrelated to whatever --annovar-protocols this run is using, so
    # any auto-downloads it triggers shouldn't land in the shared humandb dir.
    intervar_db = intervar_bin / "humandb"
    intervar_db.mkdir(parents=True, exist_ok=True)

    # InterVar's required databases are unrelated to --annovar-db and are
    # fetched lazily on first use — warn once upfront rather than let a
    # multi-GB download start silently mid-pipeline with no explanation.
    if len(list(intervar_db.glob("*.txt*"))) < 5:
        logger.warning(
            "[%s] InterVar's own reference databases aren't fully cached yet — "
            "first use may download several GB (bounded to %ds, then skips "
            "gracefully if incomplete).",
            label, _INTERVAR_TIMEOUT_S,
        )

    script = intervar_bin / "Intervar.py"
    cmd = [
        "python", str(script),
        "-b", cfg.annovar_buildver,
        "-i", str(vcf),
        "--input_type=VCF",
        "-o", str(out_prefix),
        "-d", str(intervar_db),
        "--table_annovar", str(Path(cfg.annovar_bin) / "table_annovar.pl"),
        "--convert2annovar", str(Path(cfg.annovar_bin) / "convert2annovar.pl"),
        "--annotate_variation", str(Path(cfg.annovar_bin) / "annotate_variation.pl"),
    ]

    logger.info("[%s] Running InterVar ACMG classification ...", label)
    # InterVar shells out to `annotate_variation.pl -downdb` as a *grandchild*
    # process when fetching its own missing databases. A plain subprocess.run
    # timeout only kills the direct child (Intervar.py) — the download
    # subprocess gets orphaned and keeps running (and consuming disk/network)
    # indefinitely. start_new_session + killpg ensures the whole tree dies.
    proc = subprocess.Popen(
        cmd, env=cfg.env(), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, start_new_session=True,
    )
    try:
        stdout, stderr = proc.communicate(timeout=_INTERVAR_TIMEOUT_S)
        returncode = proc.returncode
    except subprocess.TimeoutExpired:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.communicate()  # reap the killed process, avoid a zombie
        logger.warning(
            "[%s] InterVar timed out after %ds (likely downloading its own "
            "missing databases) — killed and skipping ACMG classification.",
            label, _INTERVAR_TIMEOUT_S,
        )
        return None
    if returncode != 0:
        logger.warning(
            "[%s] InterVar exited non-zero (%d) — skipping ACMG classification:\n%s",
            label, returncode, stderr[-500:],
        )
        return None

    output = Path(f"{out_prefix}.{cfg.annovar_buildver}_multianno.txt.intervar")
    if not output.exists():
        logger.warning(
            "[%s] Expected InterVar output not found at %s — skipping merge.",
            label, output,
        )
        return None
    return output


def _merge_acmg(label: str, intervar_table: Path, enriched_txt: Path) -> None:
    """Merge ACMG_classification / ACMG_evidence columns onto *enriched_txt* in place."""
    if not enriched_txt.exists():
        logger.warning(
            "[%s] Enriched annotation table not found: %s — skipping ACMG merge.",
            label, enriched_txt,
        )
        return

    # Peek at the header only (cheap) before loading the full table — the
    # InterVar output re-embeds most of ANNOVAR's own annotation columns,
    # of which only the join keys + classification column are ever used.
    # Found via audit.
    header_cols = list(pd.read_csv(intervar_table, sep="\t", nrows=0).columns)
    col = next((c for c in header_cols if _INTERVAR_COLUMN_HINT in c), None)
    if col is None:
        logger.warning(
            "[%s] Could not locate InterVar's classification column in %s — skipping merge.",
            label, intervar_table,
        )
        return
    key_cols = [c for c in ("Chr", "Start", "End", "Ref", "Alt") if c in header_cols]

    intervar = pd.read_csv(intervar_table, sep="\t", dtype=str, usecols=key_cols + [col])
    parsed = intervar[col].fillna("").str.split(n=1, expand=True)
    intervar["ACMG_classification"] = parsed[0] if 0 in parsed else ""
    intervar["ACMG_evidence"] = parsed[1] if 1 in parsed else ""

    table = pd.read_csv(enriched_txt, sep="\t", dtype=str)
    if key_cols and all(c in table.columns for c in key_cols):
        table = table.merge(
            intervar[key_cols + ["ACMG_classification", "ACMG_evidence"]],
            on=key_cols, how="left",
        )
        table.to_csv(enriched_txt, sep="\t", index=False)
        logger.log(25, "[%s] ACMG classification merged into %s", label, enriched_txt)
    else:
        logger.warning("[%s] Missing join keys — skipping ACMG merge.", label)


def run_intervar(sample: str, cfg: "Config", checkpoint: Checkpoint) -> None:
    """
    Input  : <vcf_dir>/<sample>_PASS.vcf
    Output : ACMG_classification / ACMG_evidence merged into
             <vcf_dir>/<sample>.annovar.hpo.txt
    """
    if checkpoint.done(sample, STEP):
        logger.info("[%s] ACMG classification already completed, skipping.", sample)
        return

    out_prefix = cfg.vcf_dir / f"{sample}.intervar"
    table = _run_intervar_tool(sample, cfg.vcf_dir / f"{sample}_PASS.vcf", out_prefix, cfg)
    if table is not None:
        _merge_acmg(sample, table, cfg.vcf_dir / f"{sample}.annovar.hpo.txt")
        # Only checkpoint on a genuine classification — if InterVar is
        # missing/unprovisioned/timed out, leave this un-checkpointed so a
        # later run (once InterVar/its DBs are actually available) retries
        # instead of silently skipping forever. Found via audit.
        checkpoint.mark(sample, STEP)
    else:
        logger.info(
            "[%s] ACMG classification skipped this run — will retry on next run.",
            sample,
        )


def run_intervar_cohort(samples: list[str], cfg: "Config") -> None:
    """Cohort counterpart of `run_intervar`, applied once to the cohort PASS VCF."""
    out_prefix = cfg.cohort_dir / "cohort.intervar"
    table = _run_intervar_tool(
        "cohort", cfg.cohort_dir / "cohort_PASS.vcf", out_prefix, cfg
    )
    if table is not None:
        _merge_acmg("cohort", table, cfg.cohort_dir / "cohort.annovar.hpo.txt")
