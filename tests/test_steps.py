"""
Phase 1 regression net: for every mode combination, the composed step list
must match what the V2 plan promises — most importantly that a bare
(germline, no joint-genotyping) run keeps producing the same per-sample
steps v1 always did, and cohort steps only activate when the user opts in.
"""

from exomeflow.config import Config
from exomeflow.pipeline import COHORT_STEPS, SAMPLE_STEPS, _sample_is_complete
from exomeflow.utils import Checkpoint


def _sample_names(cfg: Config) -> list[str]:
    return [s.name for s in SAMPLE_STEPS if s.applies(cfg)]


def _cohort_names(cfg: Config) -> list[str]:
    return [c.name for c in COHORT_STEPS if c.applies(cfg)]


def test_default_germline_matches_v1_step_set(cfg: Config):
    assert _sample_names(cfg) == [
        "fastp", "bwa", "sort", "flagstat", "markdup", "index", "bqsr",
        "haplotype", "filter", "annovar", "hpo", "acmg",
    ]
    assert _cohort_names(cfg) == ["multiqc"]


def test_joint_genotyping_disables_per_sample_annotation(cfg: Config):
    cfg.joint_genotyping = True
    names = _sample_names(cfg)
    assert "haplotype" in names  # still runs per-sample to produce GVCFs
    for excluded in ("filter", "annovar", "hpo", "acmg", "mutect2", "somatic_filter"):
        assert excluded not in names

    cohort_names = _cohort_names(cfg)
    assert cohort_names == [
        "joint_genotyping", "cohort_filter", "cohort_annovar",
        "cohort_hpo", "cohort_acmg", "multiqc",
    ]


def test_somatic_mode_swaps_calling_and_filtering(cfg: Config):
    cfg.mode = "somatic"
    names = _sample_names(cfg)
    assert "mutect2" in names and "haplotype" not in names
    assert "somatic_filter" in names and "filter" not in names
    # annotation/HPO/ACMG are shared across modes
    assert "annovar" in names and "hpo" in names and "acmg" in names
    assert _cohort_names(cfg) == ["multiqc"]


def test_somatic_mode_ignores_joint_genotyping(cfg: Config):
    # Cohort joint genotyping is a germline-only concept for V2 — combining
    # it with somatic mode must not silently drop per-sample somatic output.
    cfg.mode = "somatic"
    cfg.joint_genotyping = True
    names = _sample_names(cfg)
    assert "mutect2" in names and "annovar" in names
    assert _cohort_names(cfg) == ["multiqc"]


def test_cnv_step_is_independently_gated(cfg: Config):
    cfg.call_cnv = True
    names = _sample_names(cfg)
    assert "cnv" in names
    assert "filter" in names and "annovar" in names  # unaffected

    cfg.call_cnv = False
    assert "cnv" not in _sample_names(cfg)


def test_n_samples_without_joint_genotyping_yields_n_annotated_paths(cfg: Config):
    """Non-negotiable invariant from the V2 plan: no implicit cohort merge."""
    names = _sample_names(cfg)
    assert "annovar" in names
    assert not any(c.name.startswith("cohort_") for c in COHORT_STEPS if c.applies(cfg))


# --------------------------------------------------------------------------- _sample_is_complete
# Regression tests for the checkpoint-completeness redesign: completeness is
# now derived from the applicable step list, not a separate coarse flag.

def test_sample_incomplete_until_every_applicable_step_checkpointed(cfg: Config, tmp_path):
    checkpoint = Checkpoint(tmp_path / ".checkpoints")
    applicable = _sample_names(cfg)
    assert not _sample_is_complete("s1", cfg, checkpoint)

    for name in applicable[:-1]:
        checkpoint.mark("s1", name)
    assert not _sample_is_complete("s1", cfg, checkpoint)  # one step still missing

    checkpoint.mark("s1", applicable[-1])
    assert _sample_is_complete("s1", cfg, checkpoint)


def test_sample_complete_ignores_steps_that_dont_apply(cfg: Config, tmp_path):
    """cnv/mutect2/somatic_filter aren't required for a default germline run."""
    checkpoint = Checkpoint(tmp_path / ".checkpoints")
    for name in _sample_names(cfg):
        checkpoint.mark("s1", name)
    assert _sample_is_complete("s1", cfg, checkpoint)
    # Steps from other modes were never checkpointed and must not be required.
    assert not checkpoint.done("s1", "cnv")
    assert not checkpoint.done("s1", "mutect2")


def test_upgrade_scenario_new_steps_force_reprocessing(cfg: Config, tmp_path):
    """
    Regression test: found via audit — a v1-style output directory (old
    step set, no hpo/acmg) must not be permanently treated as complete once
    v2 adds new steps; only the new steps should be missing/reprocessed.
    """
    checkpoint = Checkpoint(tmp_path / ".checkpoints")
    old_v1_steps = ["fastp", "bwa", "sort", "flagstat", "markdup", "index",
                     "bqsr", "haplotype", "filter", "annovar"]
    for name in old_v1_steps:
        checkpoint.mark("s1", name)

    # hpo/acmg (new in v2) were never checkpointed for this sample.
    assert not _sample_is_complete("s1", cfg, checkpoint)

    checkpoint.mark("s1", "hpo")
    checkpoint.mark("s1", "acmg")
    assert _sample_is_complete("s1", cfg, checkpoint)
