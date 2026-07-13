# Changelog

## 2.0.0

Default-flag behavior is unchanged from v1: a bare `exomeflow run --input-dir fastq/
--output results/` still runs germline hg38 HaplotypeCaller-based calling and produces
one separate annotated file per sample, exactly as before. Everything below is opt-in.

### Added

- **Cohort joint genotyping** (`--joint-genotyping`) — GenomicsDBImport + GenotypeGVCFs
  across all samples, producing one shared cohort VCF/annotation instead of per-sample
  files. Strictly opt-in: running multiple samples without this flag still yields one
  annotated file per sample, same as v1.
- **Somatic mode** (`--mode somatic`) — tumor-only Mutect2 calling with its own
  contamination-aware filtering chain (GetPileupSummaries → CalculateContamination →
  FilterMutectCalls). Tumor-normal pairing is not yet supported.
- **Read-depth CNV calling** (`--cnv`) — GATK CollectReadCounts/DenoiseReadCounts/
  PlotDenoisedCopyRatios per sample, no panel-of-normals required.
- **GRCh37/hg19 support** (`--genome-build GRCh37`) — reference bundle, ANNOVAR
  buildver, and database downloads all follow the selected build automatically.
- **HPO + ACMG enrichment** — every annotated table is automatically joined with HPO
  gene-to-phenotype terms and ACMG/AMP classification (via InterVar), written to
  `*.annovar.hpo.txt` alongside the raw ANNOVAR output. Degrades to a skip-with-warning
  if InterVar/HPO data can't be auto-installed — it never blocks a run.
- **MultiQC cohort rollup** — `multiqc/exomeflow_report.html`, generated automatically
  at the end of every run.
- Every new tool/database (InterVar, HPO mapping, MultiQC, matplotlib, GRCh37 refs)
  follows the same auto-detect → auto-install/auto-download pattern already used for
  GATK/ANNOVAR — no new manual setup step.
- `tests/` (pytest) and GitHub Actions CI — the project had no automated tests before 2.0.

### Changed

- Pipeline steps are now composed from a pluggable registry (`exomeflow/steps.py`)
  instead of a hardcoded per-sample tuple, so germline/somatic/joint-genotyping/CNV all
  gate in and out of the same step list via `applies(cfg)` predicates.
- `matplotlib` moved to an optional extra (`pip install exomeflow[viz]`); auto-installed
  only when `--cnv` is used.

### Fixed

- Removed the duplicated dependency-check system: `exomeflow run` used to run two
  independent checkers back-to-back (`utils.check_all_requirements` and
  `setup_env.check_and_fix_dependencies`). Only the latter remains.
- Removed dead code (`utils.check_dependencies`, `utils.check_reference_files`) that was
  defined but never called.
- `matplotlib` is now declared as a real dependency instead of being silently required
  by the checker but missing from `pyproject.toml`.
- Version numbers (`pyproject.toml`, `__init__.py`, Dockerfile, Singularity def,
  CITATION.cff, .zenodo.json) are now consistent.
