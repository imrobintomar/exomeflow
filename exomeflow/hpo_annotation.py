"""
HPO gene-to-phenotype enrichment, layered onto ANNOVAR's multianno table.

Runs automatically after annotation (no opt-in flag) but degrades to a
skip-with-warning if the cached HPO mapping isn't available, rather than
failing the pipeline — annotation output without HPO columns is still useful.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from exomeflow.utils import Checkpoint

if TYPE_CHECKING:
    from exomeflow.config import Config

logger = logging.getLogger("exomeflow")

STEP = "hpo"

HPO_CACHE_DIR = Path.home() / ".exomeflow" / "hpo"
HPO_MAPPING_FILE = HPO_CACHE_DIR / "genes_to_phenotype.txt"
HPO_DOWNLOAD_URL = (
    "https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt"
)


_GENE_DIST_SUFFIX = re.compile(r"\(.*\)$")


def _split_genes(raw: object) -> list[str]:
    """
    ANNOVAR's Gene.refGene can be a single symbol, or a composite string:
    ';'-joined for variants overlapping multiple genes (e.g. "GENE1;GENE2"),
    or ','-joined with parenthesized distances for intergenic calls (e.g.
    "GENE1(12345),GENE2(6789)"). A plain merge on the raw column only ever
    matches a single canonical symbol, so multi-gene rows silently got no
    HPO terms with no signal that anything was missed. Found via audit.
    """
    if not isinstance(raw, str) or not raw:
        return []
    parts = re.split(r"[;,]", raw)
    return [_GENE_DIST_SUFFIX.sub("", p).strip() for p in parts if p.strip()]


def _load_hpo_map() -> "pd.DataFrame | None":
    if not HPO_MAPPING_FILE.exists():
        logger.warning(
            "HPO gene-to-phenotype mapping not found at %s — skipping HPO "
            "enrichment (auto-downloaded by `exomeflow run`'s dependency check).",
            HPO_MAPPING_FILE,
        )
        return None
    try:
        # Peek at the header only (cheap) to resolve the 3 columns we
        # actually need before loading the full file  this mapping file
        # has 20+ unrelated columns (disease IDs, frequencies, etc.) that
        # would otherwise be parsed and held in memory for nothing. Found
        # via audit.
        header_cols = list(pd.read_csv(HPO_MAPPING_FILE, sep="\t", comment="#", nrows=0).columns)
        gene_col = "gene_symbol" if "gene_symbol" in header_cols else header_cols[1]
        hpo_id_col = "hpo_id" if "hpo_id" in header_cols else header_cols[2]
        hpo_name_col = "hpo_name" if "hpo_name" in header_cols else header_cols[3]
        df = pd.read_csv(
            HPO_MAPPING_FILE, sep="\t", comment="#", dtype=str,
            usecols=[gene_col, hpo_id_col, hpo_name_col],
        )
    except Exception as exc:
        logger.warning("Could not parse HPO mapping file: %s", exc)
        return None

    return (
        df.groupby(gene_col)
        .agg(
            HPO_ID=(hpo_id_col, lambda s: ";".join(sorted(set(s)))),
            HPO_terms=(hpo_name_col, lambda s: ";".join(sorted(set(s)))),
        )
        .reset_index()
        .rename(columns={gene_col: "Gene.refGene"})
    )


def enrich(label: str, multianno_txt: Path, output_txt: Path) -> bool:
    """
    Merge HPO terms onto *multianno_txt* by `Gene.refGene`, write *output_txt*.

    Returns whether this call reached a genuinely complete state - used by
    the caller to decide whether to checkpoint. A missing *multianno_txt* is
    a legitimate completion (annotation gracefully skipped a 0-variant
    sample upstream - nothing to enrich, not a failure). A missing/
    unparseable HPO mapping is NOT: the table is still written without HPO
    columns (annotation output alone is useful), but this must be reported
    as incomplete so a later run retries once the mapping is available -
    found via audit: this used to always be checkpointed as done either way,
    permanently skipping enrichment for any sample that happened to run
    before the mapping was first cached.
    """
    if not multianno_txt.exists():
        logger.warning(
            "[%s] Annotated table not found: %s — skipping HPO enrichment.",
            label, multianno_txt,
        )
        return True

    hpo_map = _load_hpo_map()
    table = pd.read_csv(multianno_txt, sep="\t", dtype=str)
    if hpo_map is not None and "Gene.refGene" in table.columns:
        lookup = hpo_map.set_index("Gene.refGene")[["HPO_ID", "HPO_terms"]].to_dict("index")

        def _lookup_row(raw: object) -> "pd.Series":
            ids: set[str] = set()
            terms: set[str] = set()
            for gene in _split_genes(raw):
                hit = lookup.get(gene)
                if hit:
                    ids.update(hit["HPO_ID"].split(";"))
                    terms.update(hit["HPO_terms"].split(";"))
            return pd.Series({
                "HPO_ID": ";".join(sorted(ids)) if ids else None,
                "HPO_terms": ";".join(sorted(terms)) if terms else None,
            })

        table[["HPO_ID", "HPO_terms"]] = table["Gene.refGene"].apply(_lookup_row)
    table.to_csv(output_txt, sep="\t", index=False)

    if hpo_map is None:
        logger.info(
            "[%s] Wrote %s without HPO columns (mapping unavailable) — "
            "will retry enrichment on next run.",
            label, output_txt,
        )
        return False

    logger.success("[%s] HPO-enriched table: %s", label, output_txt)
    return True


def run_hpo_annotation(sample: str, cfg: "Config", checkpoint: Checkpoint) -> None:
    """
    Input  : <vcf_dir>/<sample>.annovar.<buildver>_multianno.txt
    Output : <vcf_dir>/<sample>.annovar.hpo.txt
    """
    if checkpoint.done(sample, STEP):
        logger.info("[%s] HPO enrichment already completed, skipping.", sample)
        return

    multianno = cfg.vcf_dir / f"{sample}.annovar.{cfg.annovar_buildver}_multianno.txt"
    complete = enrich(sample, multianno, cfg.vcf_dir / f"{sample}.annovar.hpo.txt")

    if complete:
        checkpoint.mark(sample, STEP)


def run_hpo_annotation_cohort(samples: list[str], cfg: "Config") -> bool:
    """Cohort counterpart of `run_hpo_annotation`, applied once to the cohort table."""
    multianno = cfg.cohort_dir / f"cohort.annovar.{cfg.annovar_buildver}_multianno.txt"
    return enrich("cohort", multianno, cfg.cohort_dir / "cohort.annovar.hpo.txt")
