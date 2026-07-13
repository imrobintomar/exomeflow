# Architecture

## Pluggable step registry

`exomeflow/steps.py` defines two dataclasses:

```python
@dataclass(frozen=True)
class SampleStep:
    name: str
    run: Callable[[str, Config, Checkpoint], None]
    applies: Callable[[Config], bool] = field(default=_always)

@dataclass(frozen=True)
class CohortStep:
    name: str
    run: Callable[[list[str], Config], None]
    applies: Callable[[Config], bool] = field(default=_always)
```

A `SampleStep` runs once per sample (fastp, BWA, BQSR, HaplotypeCaller,
ANNOVAR, ...). A `CohortStep` runs once after every sample's `SampleStep`s have
finished (joint genotyping, MultiQC).

`pipeline.py` builds an explicit, ordered `list[SampleStep]` and
`list[CohortStep]` from the step functions defined across the other modules
(`alignment.py`, `variant_calling.py`, `annotation.py`, `somatic.py`, `cnv.py`,
`joint_genotyping.py`, `reporting.py`, ...). There's no decorator or
auto-registration magic — each step is wired in explicitly, gated by an
`applies(cfg) -> bool` predicate that reads the `Config` flags set from CLI
options (`mode`, `genome_build`, `joint_genotyping`, `call_cnv`).

This is what makes `--mode somatic`, `--joint-genotyping`, and `--cnv`
additive instead of separate code paths: turning a flag on changes which
steps' `applies()` predicate returns `True`, not which function gets called.
For example, the per-sample filter/annotate steps and the cohort-level
filter/annotate steps are mutually exclusive via `applies()` — exactly one of
the two runs, controlled by `cfg.joint_genotyping`.

## Pipeline execution

For each sample, `run_pipeline()` iterates the ordered `SampleStep` list,
running any step whose `applies(cfg)` is `True` and skipping (via
`checkpoint.done(sample, step.name)`) any step that already completed in a
previous run. A sample is considered complete for retry purposes only when
every currently-applicable step is done — so upgrading ExomeFlow on an
existing output directory (e.g. picking up HPO/ACMG enrichment added in a
later version) correctly reprocesses just the new steps instead of skipping
the sample entirely.

After all samples finish their `SampleStep`s, the `CohortStep` list runs once,
gated the same way, under a `checkpoint.done("__cohort__", step.name)`
namespace.

## Checkpointing

Every step's completion is recorded as a `.checkpoints/<sample>.<step>.done`
marker file (process-crash-safe: a step is only marked done after it actually
finishes). Re-running the exact same command after an interruption skips every
already-completed step and resumes from where it stopped — see
[Resuming an interrupted run](usage.md#resuming-an-interrupted-run).

## Auto-provisioning

`setup_env.py` follows one pattern for every external dependency the pipeline
needs, regardless of type:

- **System binaries** (bwa, samtools, fastp, perl) — detect on `PATH`, install
  via `conda` if missing.
- **Large tools** (GATK, ANNOVAR, InterVar) — detect a bundled/known install
  location first, auto-download/clone if missing.
- **Reference data** (hg38/GRCh37 FASTA, dbSNP, Mills indels, ANNOVAR
  databases, HPO gene-to-phenotype mapping) — check for an existing file
  system-wide before downloading, so an existing install elsewhere on disk is
  reused instead of re-downloading tens of GB blind.

All of this is wired into `check_and_fix_dependencies()`, called automatically
at the start of `exomeflow run` — there's no separate manual setup step a user
has to run first, and no flag is required with no fallback.
