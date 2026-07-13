from pathlib import Path

from exomeflow.config import Config


def test_derived_dirs(cfg: Config):
    assert cfg.fastp_dir == cfg.output_dir / "filtered_fastp"
    assert cfg.map_dir == cfg.output_dir / "Mapsam"
    assert cfg.vcf_dir == cfg.output_dir / "VCF"
    assert cfg.cohort_dir == cfg.vcf_dir / "cohort"
    assert cfg.cnv_dir == cfg.output_dir / "CNV"
    assert cfg.checkpoint_dir == cfg.output_dir / ".checkpoints"


def test_defaults_match_v1_behavior(cfg: Config):
    assert cfg.mode == "germline"
    assert cfg.genome_build == "hg38"
    assert cfg.joint_genotyping is False
    assert cfg.call_cnv is False


def test_annovar_buildver_mapping(cfg: Config):
    cfg.genome_build = "hg38"
    assert cfg.annovar_buildver == "hg38"
    cfg.genome_build = "GRCh37"
    assert cfg.annovar_buildver == "hg19"


def test_grch37_construction_swaps_annovar_protocols(tmp_path: Path):
    """
    Regression test: found via audit — gnomAD v4.1 was never released for
    hg19/GRCh37, so the hg38-default protocol list must not be used
    unconditionally for a GRCh37 Config.
    """
    cfg_grch37 = Config(
        input_dir=tmp_path / "fastq", output_dir=tmp_path / "results",
        reference=tmp_path / "ref.fa", dbsnp=tmp_path / "dbsnp.vcf.gz",
        mills=tmp_path / "mills.vcf.gz", known_indels=tmp_path / "known_indels.vcf.gz",
        annovar_bin=tmp_path / "annovar", annovar_db=tmp_path / "annovar" / "humandb",
        genome_build="GRCh37",
    )
    assert "gnomad41_exome" not in cfg_grch37.annovar_protocols
    assert "gnomad211_exome" in cfg_grch37.annovar_protocols
    assert "gnomad211_genome" in cfg_grch37.annovar_protocols

    cfg_hg38 = Config(
        input_dir=tmp_path / "fastq", output_dir=tmp_path / "results",
        reference=tmp_path / "ref.fa", dbsnp=tmp_path / "dbsnp.vcf.gz",
        mills=tmp_path / "mills.vcf.gz", known_indels=tmp_path / "known_indels.vcf.gz",
        annovar_bin=tmp_path / "annovar", annovar_db=tmp_path / "annovar" / "humandb",
    )
    assert "gnomad41_exome" in cfg_hg38.annovar_protocols


def test_grch37_explicit_protocol_override_is_preserved(tmp_path: Path):
    """An explicit --annovar-protocols override must not be silently swapped."""
    custom = "refGene,avsnp150"
    cfg = Config(
        input_dir=tmp_path / "fastq", output_dir=tmp_path / "results",
        reference=tmp_path / "ref.fa", dbsnp=tmp_path / "dbsnp.vcf.gz",
        mills=tmp_path / "mills.vcf.gz", known_indels=tmp_path / "known_indels.vcf.gz",
        annovar_bin=tmp_path / "annovar", annovar_db=tmp_path / "annovar" / "humandb",
        genome_build="GRCh37", annovar_protocols=custom, annovar_operations="g,f",
    )
    assert cfg.annovar_protocols == custom


def test_setup_directories_only_creates_whats_needed(cfg: Config):
    cfg.setup_directories()
    assert cfg.vcf_dir.exists()
    assert not cfg.cohort_dir.exists()
    assert not cfg.cnv_dir.exists()

    cfg.joint_genotyping = True
    cfg.call_cnv = True
    cfg.setup_directories()
    assert cfg.cohort_dir.exists()
    assert cfg.cnv_dir.exists()
