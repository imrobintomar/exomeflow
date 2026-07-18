# Changelog

## 2.2.6

### Added

- **`--germline-resource` and `--panel-of-normals` are now auto-downloaded**
  for `--mode somatic` if not explicitly given — both are GATK's own public
  best-practices resources (no registration required, unlike ANNOVAR/OMIM),
  so `exomeflow run --mode somatic` alone is now enough to get Mutect2's
  recommended tumor-only accuracy accelerants without supplying any path.
  The germline-resource AF file is large (3 GB hg38, 14 GB GRCh37) and asks
  for confirmation before downloading (respects `--yes`); the PoN is small
  (17 MB hg38, 730 MB GRCh37) and just downloads. Passing either flag
  explicitly still overrides the auto-downloaded default. Verified against
  the live GATK bucket end-to-end — exact byte-for-byte match with the
  known file sizes.

## 2.2.5

### Added

- **`--panel-of-normals`** — pass a pre-built Panel of Normals VCF to Mutect2
  for `--mode somatic`, alongside the existing `--germline-resource`. Both
  are GATK's own documented alternative to matched tumor-normal calling
  (which ExomeFlow doesn't support yet), not a lesser substitute — a PoN
  filters recurrent sequencing artifacts that a population allele-frequency
  resource alone won't catch. Consumes an existing PoN VCF; building one
  from a set of normal samples isn't automated yet.

## 2.2.4

Root-caused by reviewing a real annotated output file: ACMG classification
had been silently failing on every single run since it was introduced —
verified fixed end-to-end against a live InterVar install.

### Fixed

- **ACMG classification silently failed on every run** — InterVar hard-
  requires `intervardb/mim2gene.txt` (a file from OMIM) to run its
  classification logic at all. Without it, InterVar prints "Error: can't
  read the OMIM file ... Please download it from http://www.omim.org/downloads"
  and never produces its output file, which the pipeline only ever saw as a
  generic "expected InterVar output not found — skipping merge" with no
  indication why. Unlike ANNOVAR, this file needs no registration — OMIM
  publishes it as a plain public download — but InterVar's own GitHub repo
  doesn't ship it, so a plain `git clone` never provisioned it. Now
  auto-downloaded during InterVar setup (and retroactively, for existing
  installs, the next time dependencies are checked).
- **InterVar was never given its own ACMG-criteria database path** — `-t`/
  `--database_intervar` (PVS1 LOF genes, PM1 domains, etc. — distinct from
  `-d`, the plain ANNOVAR annotation database) was never passed, so InterVar
  fell back to its config file's relative default path, resolved against
  the pipeline's working directory rather than InterVar's own install
  directory.

### Changed

- **InterVar now shares the main `--annovar-db` humandb directory** instead
  of maintaining a fully separate, isolated copy. Several of the databases
  InterVar's own protocol list needs (refGene, ensGene, knownGene, rmsk) are
  static, version-agnostic gene/repeat definitions that overlap with what
  the main pipeline already downloaded — reusing them avoids duplicating
  tens of GB in a second location. InterVar's other, specifically-versioned
  databases (avsnp147, dbnsfp42a, clinvar_20210501, etc. — different exact
  versions than `--annovar-protocols` uses) still download into that same
  shared directory rather than a separate one.

## 2.2.3

Three bugs found live, from a real `exomeflow setup` transcript, all in the
"do we already have this?" detection path for ANNOVAR databases.

### Fixed

- **`exomeflow setup` re-triggered the whole find-or-ask flow every time,
  ignoring a previously-saved `annovar_db`** — `run_setup()` computed its
  default database path as `annovar_bin/humandb` unconditionally, never
  checking `~/.exomeflow/config.json`'s already-saved `annovar_db` (the
  `exomeflow run` path already did this correctly; `exomeflow setup` didn't).
  A path you'd already answered interactively in an earlier run had to be
  re-entered — or, worse, re-found — on every subsequent `exomeflow setup`.
- **The same "Found N/7 databases ... missing: ..." line printed twice** —
  the default-location check tried `default_db` and `annovar_bin/humandb`
  as two separate candidates, but they're the same path whenever no
  explicit `--annovar-db` was given, so the identical check (and message)
  ran twice in a row.
- **The system-wide humandb search picked whichever hit `find` happened to
  print first, not the most complete one** — on a system with more than one
  `refGene.txt` (e.g. InterVar bundles its own small humandb subset for its
  own use), `find`'s output order isn't meaningful, so the small subset
  could be selected over the real, complete humandb. Now prefers whichever
  hit has the most `{buildver}_*.txt` files alongside it. Also: a timeout
  on that search used to discard any hits `find` had already printed before
  it ran out of time; now salvages that partial output instead.

## 2.2.2

### Fixed

- **ANNOVAR completeness checks didn't verify the paired `.idx` index
  file** — ANNOVAR's `-downdb -webfrom annovar` fetches a filter-type
  database's `.txt` data and its `.idx` index together, but a connection
  drop partway through could land the `.txt` without its `.idx`. Both
  `annovar_databases_complete()` (the pre-flight check) and
  `_step_annovar_databases()`'s download loop (which decided what to skip
  as "already present") only ever checked for the `.txt` file. A database
  missing only its index would be reported complete and never re-fetched.
  refGene is gene-based and has no `.idx` pair, so it's exempt from this
  check.

### Changed

- `exomeflow run`'s pipeline-start banner now also lists lab affiliation
  and no longer prints a GitHub link.

## 2.2.1

### Fixed

- **A partially-complete ANNOVAR humandb was reported and treated as fully
  set up** — found live: a directory with 6 of the 7 required databases
  (missing `clinvar_20240611`) was accepted as "found," and the pipeline
  proceeded straight to annotation without ever downloading the missing one.
  The completeness check only required *some* database to be present, not
  all of the required set. `_step_annovar_databases()` now uses the same
  per-file completeness check already used by the pre-flight status table,
  and a partial match falls through to the download loop for just the
  missing database(s) instead of short-circuiting as done.

## 2.2.0

Focused on the "just-installed, not-technically-sound user" path: `pip install
exomeflow` → `exomeflow setup` → `exomeflow run` with as little manual work in
between as actually achievable.

### Added

- **Auto-bootstrapped micromamba** — `bwa`/`samtools`/`fastp`/`perl` used to
  require conda or mamba already on `PATH`, with a hard stop and "go install
  Miniconda yourself" if neither was found. Now, when neither is present,
  ExomeFlow downloads a small (~7 MB) self-contained micromamba binary into
  `~/.exomeflow/conda/` and uses it to install the missing tools into an
  isolated prefix — no separate manual install step required first.
- **`exomeflow doctor`** — a new, read-only pre-flight command that reports
  every dependency's status (found/missing) plus whether it will resolve
  itself automatically or needs manual action, before any setup or download
  starts. Makes no network requests and never writes config; safe to run any
  time, including before the very first `exomeflow run`.
- `exomeflow --version` now also prints the author/institution/contact line.

### Fixed

- **ANNOVAR "not found" messaging was wrong for a real install** — it said
  "place the annovar/ folder inside the ExomeFlow directory," which is
  leftover text from a dev-checkout layout; for a `pip install`, there is no
  such directory. Corrected the message, and — matching the interactive
  prompt that already existed for the ANNOVAR *database* directory — added
  the same "where did you put it?" prompt for the ANNOVAR *scripts*
  directory, at all three places this was previously a hard, message-only
  stop (`exomeflow setup`'s pre-flight gate, `_step_bundled_tools`, and
  `check_and_fix_dependencies`'s auto-fix path). `exomeflow setup` no longer
  exits before even starting if ANNOVAR isn't auto-detected — GATK, system
  tools, and reference files still get set up regardless, since ANNOVAR is
  the one piece that can't be resolved without the user's own registered
  download.

## 2.1.3

### Added

- **`<sample>_R1.fastq.gz` / `<sample>_R2.fastq.gz` naming convention** —
  accepted alongside the existing `_1`/`_2` convention, resolved per sample,
  so a single run can mix both.

### Fixed

- **`exomeflow setup` / `exomeflow run` could "forget" a fully-configured
  ANNOVAR (and potentially GATK) install** — `detect_annovar_bin()` and
  `detect_gatk_bin()` only checked cwd-relative and hardcoded paths, never
  the path already saved in `~/.exomeflow/config.json`. Running the CLI from
  any directory other than the one it was originally set up from could
  report "ANNOVAR not found" even though it was already fully configured.
  Both detectors now check the saved config first.
- **A non-standard `DP < 10` hard filter was discarding real variants** —
  on top of GATK's own recommended germline hard-filter thresholds (QD, FS,
  SOR, MQ, MQRankSum, ReadPosRankSum — all still present, unchanged, and
  matching GATK's official guidance exactly), a site-level minimum-depth
  filter had been added that isn't part of that recommendation. Unlike the
  other thresholds, raw depth doesn't distinguish "modest but real" from
  "artifact" — QD (quality normalized by depth) already covers that — so it
  could discard genuine variants at borderline-but-real coverage (e.g. near
  capture-kit target edges). Removed; the hard-filter set is now exactly
  GATK's documented recommendation, nothing added on top.

## 2.1.2

Packaging-only fix: the README shipped as PyPI's project description in 2.1.1
rendered broken.

### Fixed

- **Badges table rendered as literal `|---|---|` text on PyPI** — the table's
  delimiter row had no header row above it, which PyPI's stricter README
  renderer requires (GitHub tolerated it). Added the missing header row.
- **Logo didn't render on PyPI** — the `<img>` used a relative `src` path,
  which resolves against the GitHub repo on GitHub's renderer but has no
  equivalent base on PyPI. Switched to the `raw.githubusercontent.com`
  absolute URL.

No code changes — `exomeflow` behavior is identical to 2.1.1.

## 2.1.1

Bug-fix release following a full source audit of 2.0.0. The most important fix:
**`--joint-genotyping` was completely broken in 2.0.0** — every cohort run crashed
after all the expensive GATK work had already succeeded. If you use cohort mode,
upgrade before relying on it.

### Fixed

- **`--joint-genotyping` crashed every run** — `count_variants()` opened every VCF as
  plain UTF-8 text with no gzip handling, but the cohort path feeds it GATK's
  bgzipped `cohort.vcf.gz`. The resulting `UnicodeDecodeError` wasn't a
  `PipelineStepError`, so it went uncaught and crashed the whole process after every
  GATK step had already succeeded.
- **Switching `--genome-build` on an already-configured install silently kept the old
  build's reference/ANNOVAR files** — the check was pure file-existence with no
  comparison against the previously-saved build, so BWA/GATK could run against hg38
  while ANNOVAR/InterVar annotated as if the variants were hg19, with no error.
- **GRCh37 annotation was broken even with correct references** — the default ANNOVAR
  protocol list (`gnomad41_exome`/`gnomad41_genome`) is hg38-only; gnomAD v4.1 was
  never released for hg19/GRCh37. GRCh37 runs now correctly use `gnomad211_exome`/
  `gnomad211_genome` instead. `--genome-build` also now correctly persists to
  `~/.exomeflow/config.json` on every path, not just when something needed fixing.
- **Per-tool minimum-version enforcement was silently dropped** when the two
  dependency checkers were consolidated in 2.0.0 — only presence was checked, not
  actual version. Restored: an outdated tool is now detected and conda-upgraded
  automatically, with a hard failure if the upgrade doesn't clear the minimum.
- **ANNOVAR database completeness narrowed to "does the directory exist"** — a
  database deleted after initial setup (or simply missing for the requested build)
  went undetected until `table_annovar.pl` failed hours into a real run. Restored to
  a per-database-file check, and made build-aware (hg38 vs. GRCh37 need different
  database sets).
- **A sample could get permanently stuck unable to retry a skipped step** (e.g. ACMG
  classification skipped because InterVar wasn't provisioned yet) — the coarse
  per-sample "COMPLETE" checkpoint flag was set regardless, so even after fixing the
  underlying cause, the sample was never reprocessed. Checkpoint completeness (for
  retry purposes) is now derived from the actual applicable step list instead of a
  separate flag — this also means upgrading an existing v1/v2.0.0 output directory
  correctly reprocesses just the new steps (hpo/acmg) instead of skipping the sample
  entirely. (A related regression introduced while fixing this — a gracefully-skipped
  optional step incorrectly excluded the sample from the MultiQC cohort rollup even on
  an otherwise fully successful run — was caught by live smoke testing before release
  and fixed in the same pass.)
- **Several subprocess calls could hang forever** on a stalled connection, including
  the exact script (`annotate_variation.pl -downdb`) whose orphaned-child-process
  behavior had already been fixed for InterVar in 2.0.0 but not applied to its other
  call sites. All now have bounded, process-group-killing timeouts; large file
  downloads (wget/curl) use stall-detection instead of a fixed cutoff so legitimate
  multi-hour transfers aren't killed.
- **The Python version floor check (<3.9) was dropped entirely** with no replacement
  — restored with a clear upfront error instead of a confusing failure deep inside a
  3.9+-only module.
- `shutil.disk_usage()` crash on system-resource auto-detection when `--output`
  doesn't exist yet (the common case — it's documented as "created if absent").
- `--annovar-protocols`/`--annovar-operations` CLI overrides added — a user with a
  pre-existing, differently-versioned ANNOVAR install now has a way to work around a
  version mismatch without editing source.
- Fail-fast validation moved before the (potentially multi-hour) setup wizard: a bad
  `--input-dir`, or `--joint-genotyping`/`--cnv` without `--intervals`, now fails in
  under a second instead of after the wizard completes.
- ANNOVAR crashing on a completely empty PASS VCF (e.g. catastrophic sequencing
  failure, or an intervals BED with near-zero coverage overlap) — now skips annotation
  cleanly with a warning instead of a cryptic Perl error.
- System resource auto-detection (`--threads`/`--java-opts` sizing from actual
  CPU/RAM) and a system-wide ANNOVAR humandb lookup (reuses an existing installation
  elsewhere on disk instead of re-downloading ~90GB blind).

### Changed

- Deduplicated repeated logic found during the audit: the hg38/hg19 ANNOVAR-buildver
  mapping (previously computed independently in three places), the genome-build CLI
  validation (duplicated in `run`/`setup`), reference-file alternate-name lookup
  (duplicated as a nested function), and the `--intervals`-presence check (reimplemented
  in `cli.py` instead of reusing `Config.has_intervals`).
- ANNOVAR protocol/operation pairing is now validated in `Config` itself (previously
  only in `cli.py`, so any other caller building a `Config` directly bypassed it).
- HPO/ACMG pandas reads now load only the columns actually used instead of the full
  table.

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
