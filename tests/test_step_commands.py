"""
Command-argv construction tests. Nothing here executes a real subprocess —
`run_cmd` (or, where a module calls subprocess directly, subprocess.run /
Popen) is monkeypatched to a recorder so we can assert on the built argv
without needing bwa/GATK/ANNOVAR/InterVar actually installed.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from exomeflow.config import Config
from exomeflow.utils import Checkpoint


def _recorder(calls: list):
    def _fake(cmd, **kwargs):
        calls.append(cmd)
        return MagicMock(returncode=0)
    return _fake


@pytest.fixture
def checkpoint(tmp_path: Path) -> Checkpoint:
    return Checkpoint(tmp_path / ".checkpoints")


def _touch_fastq(cfg: Config, sample: str):
    cfg.input_dir.mkdir(parents=True, exist_ok=True)
    (cfg.input_dir / f"{sample}_1.fastq.gz").touch()
    (cfg.input_dir / f"{sample}_2.fastq.gz").touch()


# --------------------------------------------------------------------------- fastp

def test_fastp_command(cfg: Config, checkpoint: Checkpoint, monkeypatch):
    import exomeflow.fastp as mod
    _touch_fastq(cfg, "s1")
    calls: list = []
    monkeypatch.setattr(mod, "run_cmd", _recorder(calls))
    mod.run_fastp("s1", cfg, checkpoint)
    assert calls[0][0] == "fastp"
    assert "--length_required" in calls[0] and "50" in calls[0]
    assert checkpoint.done("s1", "fastp")


# --------------------------------------------------------------------------- bwa

def test_bwa_pipes_into_samtools(cfg: Config, checkpoint: Checkpoint, monkeypatch):
    import exomeflow.alignment as mod

    popen_calls: list = []

    class FakeProc:
        def __init__(self, cmd, **kwargs):
            popen_calls.append(cmd)
            self.stdout = MagicMock()
            self.returncode = 0

        def wait(self):
            return 0

    monkeypatch.setattr(mod.subprocess, "Popen", FakeProc)
    mod.run_bwa_mem("s1", cfg, checkpoint)
    assert popen_calls[0][:2] == ["bwa", "mem"]
    assert popen_calls[1][:2] == ["samtools", "view"]
    assert checkpoint.done("s1", "bwa")


# --------------------------------------------------------------------------- bam_processing

def test_sort_and_markdup_and_index_commands(cfg: Config, checkpoint: Checkpoint, monkeypatch):
    import exomeflow.bam_processing as mod
    calls: list = []
    monkeypatch.setattr(mod, "run_cmd", _recorder(calls))

    mod.sort_bam("s1", cfg, checkpoint)
    mod.mark_duplicates("s1", cfg, checkpoint)
    mod.build_bam_index("s1", cfg, checkpoint)

    assert calls[0][:2] == ["gatk", "SortSam"]
    assert calls[1][:2] == ["gatk", "MarkDuplicates"]
    assert calls[2][:2] == ["gatk", "BuildBamIndex"]


def test_flagstat_writes_stdout_to_file(cfg: Config, checkpoint: Checkpoint, monkeypatch):
    import exomeflow.bam_processing as mod
    cfg.setup_directories()

    def fake_run(cmd, **kwargs):
        assert cmd == ["samtools", "flagstat", str(cfg.map_dir / "s1_sorted.bam")]
        return MagicMock(returncode=0, stdout="42 mapped\n", stderr="")

    monkeypatch.setattr(mod.subprocess, "run", fake_run)
    mod.generate_flagstat("s1", cfg, checkpoint)
    assert (cfg.map_dir / "s1_flagstat.txt").read_text() == "42 mapped\n"


# --------------------------------------------------------------------------- bqsr

def test_bqsr_known_sites(cfg: Config, checkpoint: Checkpoint, monkeypatch):
    import exomeflow.recalibration as mod
    calls: list = []
    monkeypatch.setattr(mod, "run_cmd", _recorder(calls))
    mod.run_base_recalibration("s1", cfg, checkpoint)
    recal_cmd = calls[0]
    assert recal_cmd[:2] == ["gatk", "BaseRecalibrator"]
    assert recal_cmd.count("--known-sites") == 3
    assert calls[1][:2] == ["gatk", "ApplyBQSR"]


# --------------------------------------------------------------------------- variant_calling (GVCF branch)

def test_haplotype_caller_default_vs_gvcf(cfg: Config, checkpoint: Checkpoint, monkeypatch):
    import exomeflow.variant_calling as mod
    calls: list = []
    monkeypatch.setattr(mod, "run_cmd", _recorder(calls))

    mod.run_haplotype_caller("s1", cfg, checkpoint)
    assert "-ERC" not in calls[0]
    assert str(cfg.vcf_dir / "s1.vcf") in calls[0]

    checkpoint2 = Checkpoint(cfg.checkpoint_dir.parent / ".checkpoints2")
    cfg.joint_genotyping = True
    mod.run_haplotype_caller("s1", cfg, checkpoint2)
    assert "-ERC" in calls[1] and "GVCF" in calls[1]
    assert str(cfg.vcf_dir / "s1.g.vcf.gz") in calls[1]


# --------------------------------------------------------------------------- filtering (shared hard_filter)

def test_hard_filter_shared_by_sample_and_cohort_paths(cfg: Config, monkeypatch):
    import exomeflow.filtering as mod
    calls: list = []
    monkeypatch.setattr(mod, "run_cmd", _recorder(calls))
    cfg.setup_directories()

    mod.hard_filter(
        label="s1", input_vcf=cfg.vcf_dir / "s1.vcf",
        output_pass_vcf=cfg.vcf_dir / "s1_PASS.vcf", workdir=cfg.vcf_dir,
        reference=cfg.reference, env=cfg.env(),
    )
    sample_steps = [c[1] for c in calls]  # gatk subcommand names

    calls.clear()
    cfg.joint_genotyping = True
    cfg.setup_directories()
    mod.hard_filter(
        label="cohort", input_vcf=cfg.cohort_dir / "cohort.vcf.gz",
        output_pass_vcf=cfg.cohort_dir / "cohort_PASS.vcf", workdir=cfg.cohort_dir,
        reference=cfg.reference, env=cfg.env(),
    )
    cohort_steps = [c[1] for c in calls]

    assert sample_steps == cohort_steps == [
        "SelectVariants", "SelectVariants", "VariantFiltration",
        "VariantFiltration", "MergeVcfs", "SelectVariants",
    ]


# --------------------------------------------------------------------------- annotation (buildver-aware)

_FAKE_VCF_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _write_fake_pass_vcf(cfg: Config, sample: str, n_variants: int = 1) -> None:
    cfg.vcf_dir.mkdir(parents=True, exist_ok=True)
    lines = _FAKE_VCF_HEADER
    for i in range(n_variants):
        lines += f"chr1\t{100 + i}\t.\tA\tG\t50\tPASS\t.\n"
    (cfg.vcf_dir / f"{sample}_PASS.vcf").write_text(lines)


def test_annotation_uses_annovar_buildver(cfg: Config, checkpoint: Checkpoint, monkeypatch):
    import exomeflow.annotation as mod
    calls: list = []
    monkeypatch.setattr(mod, "run_cmd", _recorder(calls))
    _write_fake_pass_vcf(cfg, "s1")

    cfg.genome_build = "GRCh37"
    mod.run_annovar_annotation("s1", cfg, checkpoint)
    cmd = calls[0]
    assert "--buildver" in cmd
    assert cmd[cmd.index("--buildver") + 1] == "hg19"


def test_annotation_skips_cleanly_on_zero_variants(cfg: Config, checkpoint: Checkpoint, monkeypatch):
    """
    Regression test: found live during smoke testing — table_annovar.pl
    crashes ungracefully on a completely empty VCF (a real sample can hit
    this on catastrophic sequencing failure, not just a smoke-test artifact).
    """
    import exomeflow.annotation as mod
    calls: list = []
    monkeypatch.setattr(mod, "run_cmd", _recorder(calls))
    _write_fake_pass_vcf(cfg, "s1", n_variants=0)

    mod.run_annovar_annotation("s1", cfg, checkpoint)  # must not raise
    assert calls == []  # never invoked table_annovar.pl on an empty VCF
    assert checkpoint.done("s1", "annovar")  # still marked complete, not failed


# --------------------------------------------------------------------------- somatic

def test_mutect2_and_somatic_filter_commands(cfg: Config, checkpoint: Checkpoint, monkeypatch):
    import exomeflow.somatic as mod
    calls: list = []
    monkeypatch.setattr(mod, "run_cmd", _recorder(calls))

    mod.run_mutect2("s1", cfg, checkpoint)
    assert calls[0][:2] == ["gatk", "Mutect2"]
    assert "-tumor" in calls[0] and "s1" in calls[0]

    calls.clear()
    mod.run_somatic_filtration("s1", cfg, checkpoint)
    assert calls[-2][:2] == ["gatk", "FilterMutectCalls"]
    assert calls[-1][:2] == ["gatk", "SelectVariants"]
    # no --germline-resource supplied -> no pileup/contamination calls
    assert not any(c[:2] == ["gatk", "GetPileupSummaries"] for c in calls)


def test_mutect2_uses_germline_resource_when_supplied(cfg: Config, checkpoint: Checkpoint, monkeypatch):
    import exomeflow.somatic as mod
    calls: list = []
    monkeypatch.setattr(mod, "run_cmd", _recorder(calls))
    cfg.germline_resource = Path("/refs/af-only-gnomad.vcf.gz")

    mod.run_mutect2("s1", cfg, checkpoint)
    assert "--germline-resource" in calls[0]

    calls.clear()
    mod.run_somatic_filtration("s1", cfg, checkpoint)
    assert any(c[:2] == ["gatk", "GetPileupSummaries"] for c in calls)
    assert any(c[:2] == ["gatk", "CalculateContamination"] for c in calls)


def test_mutect2_uses_panel_of_normals_when_supplied(cfg: Config, checkpoint: Checkpoint, monkeypatch):
    import exomeflow.somatic as mod
    calls: list = []
    monkeypatch.setattr(mod, "run_cmd", _recorder(calls))
    cfg.panel_of_normals = Path("/refs/pon.vcf.gz")

    mod.run_mutect2("s1", cfg, checkpoint)
    assert "--panel-of-normals" in calls[0]
    assert calls[0][calls[0].index("--panel-of-normals") + 1] == "/refs/pon.vcf.gz"


def test_mutect2_omits_panel_of_normals_when_not_supplied(cfg: Config, checkpoint: Checkpoint, monkeypatch):
    import exomeflow.somatic as mod
    calls: list = []
    monkeypatch.setattr(mod, "run_cmd", _recorder(calls))
    assert cfg.panel_of_normals is None

    mod.run_mutect2("s1", cfg, checkpoint)
    assert "--panel-of-normals" not in calls[0]


# --------------------------------------------------------------------------- CNV

def test_cnv_requires_intervals(cfg: Config, checkpoint: Checkpoint):
    from exomeflow.cnv import run_cnv_calling
    from exomeflow.utils import PipelineStepError

    with pytest.raises(PipelineStepError):
        run_cnv_calling("s1", cfg, checkpoint)


def test_cnv_command_chain(cfg: Config, checkpoint: Checkpoint, monkeypatch, tmp_path: Path):
    import exomeflow.cnv as mod
    calls: list = []
    monkeypatch.setattr(mod, "run_cmd", _recorder(calls))
    cfg.intervals = tmp_path / "targets.bed"
    cfg.intervals.touch()

    mod.run_cnv_calling("s1", cfg, checkpoint)
    subcommands = [c[1] for c in calls]
    assert subcommands == ["CollectReadCounts", "DenoiseReadCounts", "PlotDenoisedCopyRatios"]


# --------------------------------------------------------------------------- joint genotyping

def test_joint_genotyping_requires_intervals(cfg: Config):
    from exomeflow.joint_genotyping import run_joint_genotyping
    from exomeflow.utils import PipelineStepError

    with pytest.raises(PipelineStepError):
        run_joint_genotyping(["s1", "s2"], cfg)


def test_joint_genotyping_command_chain(cfg: Config, monkeypatch, tmp_path: Path):
    import exomeflow.joint_genotyping as mod
    calls: list = []
    monkeypatch.setattr(mod, "run_cmd", _recorder(calls))
    cfg.intervals = tmp_path / "targets.bed"
    cfg.intervals.touch()

    mod.run_joint_genotyping(["s1", "s2"], cfg)
    assert calls[0][:2] == ["gatk", "GenomicsDBImport"]
    assert calls[0].count("-V") == 2
    assert calls[1][:2] == ["gatk", "GenotypeGVCFs"]


# --------------------------------------------------------------------------- multiqc (graceful skip)

def test_multiqc_skips_when_not_on_path(cfg: Config, monkeypatch):
    import exomeflow.reporting as mod
    monkeypatch.setattr(mod.shutil, "which", lambda name: None)
    called = []
    monkeypatch.setattr(mod.subprocess, "run", lambda *a, **k: called.append(a))
    mod.run_multiqc(["s1"], cfg)
    assert called == []  # never invoked subprocess when multiqc isn't installed
