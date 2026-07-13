"""
HPO gene-to-phenotype enrichment, layered onto ANNOVAR's multianno table.

Runs automatically after annotation (no opt-in flag) but degrades to a
skip-with-warning if the cached HPO mapping isn't available, rather than
failing the pipeline — annotation output without HPO columns is still useful.
"""

from __future__ import annotations

import logging
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
        # actually need before loading the full file — this mapping file
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


def enrich(label: str, multianno_txt: Path, output_txt: Path) -> None:
    """Merge HPO terms onto *multianno_txt* by `Gene.refGene`, write *output_txt*."""
    if not multianno_txt.exists():
        logger.warning(
            "[%s] Annotated table not found: %s — skipping HPO enrichment.",
            label, multianno_txt,
        )
        return

    hpo_map = _load_hpo_map()
    table = pd.read_csv(multianno_txt, sep="\t", dtype=str)
    if hpo_map is not None and "Gene.refGene" in table.columns:
        table = table.merge(hpo_map, on="Gene.refGene", how="left")
    table.to_csv(output_txt, sep="\t", index=False)
    logger.log(25, "[%s] HPO-enriched table: %s", label, output_txt)


def run_hpo_annotation(sample: str, cfg: "Config", checkpoint: Checkpoint) -> None:
    """
    Input  : <vcf_dir>/<sample>.annovar.<buildver>_multianno.txt
    Output : <vcf_dir>/<sample>.annovar.hpo.txt
    """
    if checkpoint.done(sample, STEP):
        logger.info("[%s] HPO enrichment already completed, skipping.", sample)
        return

    multianno = cfg.vcf_dir / f"{sample}.annovar.{cfg.annovar_buildver}_multianno.txt"
    enrich(sample, multianno, cfg.vcf_dir / f"{sample}.annovar.hpo.txt")

    checkpoint.mark(sample, STEP)


def run_hpo_annotation_cohort(samples: list[str], cfg: "Config") -> None:
    """Cohort counterpart of `run_hpo_annotation`, applied once to the cohort table."""
    multianno = cfg.cohort_dir / f"cohort.annovar.{cfg.annovar_buildver}_multianno.txt"
    enrich("cohort", multianno, cfg.cohort_dir / "cohort.annovar.hpo.txt")
