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


# InterVar's config.ini pins its own specifically-versioned DB set
# (avsnp147, dbnsfp42a, clinvar_20210501, 1000g2015aug, esp6500siv2_all,
# gnomad_genome, dbscsnv11, ...) that mostly doesn't overlap by exact
# filename with whatever --annovar-protocols downloaded, even though both
# now share the same humandb directory (2.2.6) — so a first-ever ACMG
# classification on a machine can still trigger InterVar auto-downloading
# tens of GB of its own databases (dbnsfp42a alone is ~48GB). 1800s (30 min)
# was nowhere near enough for that and left ACMG stuck re-timing-out on
# every run — found live. Matched to the 10800s already used for a single
# ANNOVAR -downdb call, doubled since InterVar can need several such files
# in one pass, not just one.
_INTERVAR_TIMEOUT_S = 21600


def _intervar_required_db_names(intervar_bin: Path) -> list[str]:
    """
    Parse the `database_names` line out of InterVar's own config.ini,
    rather than hardcoding its DB set here — that list is InterVar's own
    to define and can shift across InterVar versions/installs.
    """
    try:
        text = (intervar_bin / "config.ini").read_text()
    except OSError:
        return []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("database_names") and "=" in line:
            names = line.split("=", 1)[1].split()
            # "1000g2015aug" isn't a literal filename - it expands to
            # per-population site files.
            expanded: list[str] = []
            for n in names:
                if n == "1000g2015aug":
                    expanded += [f"{pop}.sites.2015_08" for pop in
                                 ("AFR", "AMR", "EAS", "EUR", "SAS", "ALL")]
                else:
                    expanded.append(n)
            return expanded
    return []


def _resolve_intervar_db(intervar_bin: Path, cfg: "Config") -> Path:
    """
    Pick whichever candidate directory already has more of InterVar's own
    required database set, instead of always preferring the shared
    --annovar-db. Found live: the shared humandb tracks the main pipeline's
    newer database versions (e.g. avsnp150, dbnsfp47a), but InterVar's
    config.ini hardcodes older, specific versions (avsnp147, dbnsfp42a,
    clinvar_20210501, ...) that don't exist under those newer names. A
    machine with a legacy standalone InterVar install can have its own
    humandb/ already fully populated with exactly what InterVar wants —
    unconditionally preferring the shared dir meant that complete, already-
    downloaded set was ignored, and InterVar instead attempted to
    re-download a ~48GB dbnsfp42a from ANNOVAR's own server (observed live
    to run at ~500 bytes/sec — a multi-year download at that rate).
    """
    own_db = intervar_bin / "humandb"
    shared_db = Path(cfg.annovar_db) if cfg.annovar_db else None
    required = _intervar_required_db_names(intervar_bin)
    buildver = cfg.annovar_buildver

    def _present_count(d: Path | None) -> int:
        if not d or not d.is_dir():
            return 0
        return sum(1 for name in required if (d / f"{buildver}_{name}.txt").exists())

    own_count = _present_count(own_db)
    shared_count = _present_count(shared_db)
    if own_count >= shared_count and own_count > 0:
        return own_db
    return shared_db if shared_db else own_db


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

    # InterVar's own protocol list (refGene, esp6500siv2_all, 1000g2015aug,
    # avsnp147, dbnsfp42a, clinvar_20210501, gnomad_genome, dbscsnv11, rmsk,
    # ensGene, knownGene) is fixed by its own config.ini and mostly doesn't
    # match --annovar-protocols' versions — but several of these (refGene,
    # ensGene, knownGene, rmsk) are static, version-agnostic gene/repeat
    # definitions that overlap with what the main pipeline already
    # downloaded. Point -d at the same shared humandb instead of a second,
    # fully isolated copy: whatever already overlaps is reused for free, and
    # whatever InterVar still needs downloads into the same shared pool
    # instead of duplicating tens of GB in a second location.
    intervar_db = _resolve_intervar_db(intervar_bin, cfg)
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
        # InterVar's own ACMG-criteria database (PVS1 LOF genes, PM1 domains,
        # OMIM mim2gene.txt, etc.) — distinct from -d, which is only the
        # ANNOVAR annotation database. Omitting this flag left InterVar
        # falling back to its config.ini's relative default ("intervardb"),
        # resolved against whatever the pipeline's cwd happened to be rather
        # than InterVar's own install directory. Found via a live run.
        "-t", str(intervar_bin / "intervardb"),
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
        # Found live: InterVar can exit 0 with no output at all if its
        # *internal* table_annovar.pl call fails (e.g. can't produce
        # `{out_prefix}.{buildver}_multianno.txt` for InterVar's own glob
        # to find) — Intervar.py doesn't propagate that as a non-zero
        # process exit. Previously this branch discarded stdout/stderr
        # entirely, so a returncode-0-but-empty-handed run was completely
        # silent about the real cause. Now surfaces the tail of both so the
        # next run's log actually explains what happened.
        logger.warning(
            "[%s] Expected InterVar output not found at %s — skipping merge. "
            "InterVar's own output (last 1500 chars of stdout+stderr):\n%s",
            label, output, (stdout[-1500:] + stderr[-1500:]).strip() or "(no output captured)",
        )
        return None
    return output


def _merge_acmg(label: str, intervar_table: Path, enriched_txt: Path) -> bool:
    """
    Merge ACMG_classification / ACMG_evidence columns onto *enriched_txt* in
    place. Returns whether columns were actually written — used by the
    caller to decide whether to checkpoint. Every early-return path here
    used to be silently treated as success by the caller (InterVar produced
    *a* file, so `run_intervar` checkpointed the step) even when this
    function did nothing — e.g. InterVar's classification column not found
    via the `_INTERVAR_COLUMN_HINT` substring match (a real risk: this
    module's own header notes InterVar's exact output columns were never
    hands-on verified). Found via audit — same "checkpoint on a graceful
    skip" bug class already fixed for HPO/MultiQC, missed here.
    """
    if not enriched_txt.exists():
        logger.warning(
            "[%s] Enriched annotation table not found: %s — skipping ACMG merge.",
            label, enriched_txt,
        )
        return False

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
        return False

    # InterVar's own first column is literally "#Chr" (leading '#'), not
    # "Chr" like the main enriched table and every other column here.
    # Found live: matching on the bare "Chr" name silently dropped
    # chromosome out of the merge join key entirely.
    chr_col = next((c for c in header_cols if c.lstrip("#") == "Chr"), None)
    key_cols = [c for c in ("Start", "End", "Ref", "Alt") if c in header_cols]
    read_cols = key_cols + [col] + ([chr_col] if chr_col else [])

    intervar = pd.read_csv(intervar_table, sep="\t", dtype=str, usecols=read_cols)
    if chr_col:
        intervar = intervar.rename(columns={chr_col: "Chr"})
        key_cols = ["Chr"] + key_cols

    # Found live: InterVar's raw value is
    # " InterVar: <classification words> PVS1=... PS=[...] ...". Splitting
    # on the first whitespace (the old approach) grabbed the literal label
    # "InterVar:" as the classification and left the real (often
    # multi-word, e.g. "Likely benign") text stuck inside "evidence"
    # instead — every row's ACMG_classification was the same useless
    # literal string "InterVar:".
    extracted = intervar[col].fillna("").str.extract(
        r"InterVar:\s*(?P<classification>.*?)\s*(?P<evidence>PVS1=.*)"
    )
    intervar["ACMG_classification"] = extracted["classification"].str.strip()
    intervar["ACMG_evidence"] = extracted["evidence"].str.strip()

    # Found live: InterVar emits one row per overlapping transcript/gene
    # for the same variant (standard ANNOVAR gene-based annotation
    # behavior for multi-isoform genes) — its classification is identical
    # across those duplicate rows for a given variant (verified against a
    # real run), but merging without deduplicating first multiplies rows
    # in the output 1:N instead of 1:1 (observed live: 937 real variants
    # became 1261 output rows).
    if key_cols:
        intervar = intervar.drop_duplicates(subset=key_cols, keep="first")

    table = pd.read_csv(enriched_txt, sep="\t", dtype=str)
    if key_cols and all(c in table.columns for c in key_cols):
        left = table.copy()
        right = intervar[key_cols + ["ACMG_classification", "ACMG_evidence"]].copy()
        merge_on = list(key_cols)
        if "Chr" in key_cols:
            # Found live: InterVar's own Chr values are bare numbers
            # ("1", "2", ..., "X") while the main annotated table uses
            # UCSC-style "chr"-prefixed values ("chr1", "chr2", ...).
            # Merging on the literal value (as introduced by the #Chr fix
            # above) matched almost nothing — of 937 real variants, only
            # 4 matched, and those were HLA allele contigs that happen to
            # already lack a "chr" prefix in both tools' output. Normalize
            # both sides for the join only; the displayed Chr value in the
            # output stays exactly as the main table had it.
            left["_chr_norm"] = left["Chr"].str.replace(r"(?i)^chr", "", regex=True)
            right["_chr_norm"] = right["Chr"].str.replace(r"(?i)^chr", "", regex=True)
            right = right.drop(columns=["Chr"])
            merge_on = ["_chr_norm" if c == "Chr" else c for c in merge_on]
        merged = left.merge(right, on=merge_on, how="left")
        if "_chr_norm" in merged.columns:
            merged = merged.drop(columns=["_chr_norm"])
        merged.to_csv(enriched_txt, sep="\t", index=False)
        logger.success("[%s] ACMG classification merged into %s", label, enriched_txt)
        return True
    logger.warning("[%s] Missing join keys — skipping ACMG merge.", label)
    return False


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
    merged = table is not None and _merge_acmg(sample, table, cfg.vcf_dir / f"{sample}.annovar.hpo.txt")
    # Only checkpoint on a genuine classification — if InterVar is
    # missing/unprovisioned/timed out, or ran but the merge itself couldn't
    # locate/attach ACMG columns, leave this un-checkpointed so a later run
    # retries instead of silently skipping forever. Found via audit: the
    # merge outcome wasn't checked at all — only whether InterVar produced
    # *a* file — so a merge failure was checkpointed as done anyway.
    if merged:
        checkpoint.mark(sample, STEP)
    else:
        logger.info(
            "[%s] ACMG classification skipped this run — will retry on next run.",
            sample,
        )


def run_intervar_cohort(samples: list[str], cfg: "Config") -> bool:
    """
    Cohort counterpart of `run_intervar`, applied once to the cohort PASS VCF.

    Returns whether ACMG columns were actually merged — pipeline.py's
    cohort loop only checkpoints on a non-False return, so InterVar being
    missing/unprovisioned (or the merge failing) correctly leaves this
    retryable instead of checkpointed done. Found via audit: this used to
    always implicitly return None, which pipeline.py treats as success.
    """
    out_prefix = cfg.cohort_dir / "cohort.intervar"
    table = _run_intervar_tool(
        "cohort", cfg.cohort_dir / "cohort_PASS.vcf", out_prefix, cfg
    )
    if table is None:
        return False
    return _merge_acmg("cohort", table, cfg.cohort_dir / "cohort.annovar.hpo.txt")
