from pathlib import Path

import exomeflow.hpo_annotation as mod


def test_enrich_merges_hpo_terms_by_gene(tmp_path: Path, monkeypatch):
    mapping = tmp_path / "genes_to_phenotype.txt"
    mapping.write_text(
        "ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\n"
        "1\tBRCA1\tHP:0003002\tBreast carcinoma\n"
        "1\tBRCA1\tHP:0100615\tOvarian neoplasm\n"
        "2\tTP53\tHP:0002664\tNeoplasm\n"
    )
    monkeypatch.setattr(mod, "HPO_MAPPING_FILE", mapping)

    multianno = tmp_path / "s1.annovar.hg38_multianno.txt"
    multianno.write_text("Chr\tStart\tGene.refGene\n1\t100\tBRCA1\n2\t200\tUNKNOWN_GENE\n")
    output = tmp_path / "s1.annovar.hpo.txt"

    complete = mod.enrich("s1", multianno, output)
    assert complete is True

    lines = output.read_text().splitlines()
    header = lines[0].split("\t")
    assert "HPO_ID" in header and "HPO_terms" in header

    brca1_row = next(row for row in lines[1:] if "BRCA1" in row)
    assert "HP:0003002" in brca1_row and "HP:0100615" in brca1_row

    unknown_row = next(row for row in lines[1:] if "UNKNOWN_GENE" in row)
    # pandas writes NaN for unmatched merges
    assert unknown_row.split("\t")[-1] in ("", "nan")


def test_enrich_skips_gracefully_when_mapping_missing(tmp_path: Path, monkeypatch):
    """
    Regression test: the table must still be written (annotation output
    without HPO columns is still useful) but this must NOT be reported as
    complete — found via audit, this used to be indistinguishable from a
    real success, so a sample processed before the HPO mapping was first
    cached could never retry enrichment once the mapping became available.
    """
    monkeypatch.setattr(mod, "HPO_MAPPING_FILE", tmp_path / "does_not_exist.txt")

    multianno = tmp_path / "s1.annovar.hg38_multianno.txt"
    multianno.write_text("Chr\tStart\tGene.refGene\n1\t100\tBRCA1\n")
    output = tmp_path / "s1.annovar.hpo.txt"

    complete = mod.enrich("s1", multianno, output)  # must not raise
    assert complete is False
    assert output.exists()
    assert "HPO_ID" not in output.read_text().splitlines()[0]


def test_enrich_skips_when_multianno_missing(tmp_path: Path):
    """A 0-variant sample where annotation was itself gracefully skipped
    upstream is a legitimate completion, not a failure to retry."""
    output = tmp_path / "out.txt"
    complete = mod.enrich("s1", tmp_path / "missing.txt", output)
    assert complete is True
    assert not output.exists()


def test_run_hpo_annotation_does_not_checkpoint_when_mapping_missing(tmp_path: Path, monkeypatch):
    from exomeflow.config import Config
    from exomeflow.utils import Checkpoint

    monkeypatch.setattr(mod, "HPO_MAPPING_FILE", tmp_path / "does_not_exist.txt")

    cfg = Config(
        input_dir=tmp_path / "fastq",
        output_dir=tmp_path / "results",
        reference=tmp_path / "ref.fa",
        dbsnp=tmp_path / "dbsnp.vcf.gz",
        mills=tmp_path / "mills.vcf.gz",
        known_indels=tmp_path / "known_indels.vcf.gz",
        annovar_bin=tmp_path / "annovar",
        annovar_db=tmp_path / "annovar" / "humandb",
    )
    cfg.vcf_dir.mkdir(parents=True, exist_ok=True)
    multianno = cfg.vcf_dir / f"s1.annovar.{cfg.annovar_buildver}_multianno.txt"
    multianno.write_text("Chr\tStart\tGene.refGene\n1\t100\tBRCA1\n")

    checkpoint = Checkpoint(cfg.checkpoint_dir)
    mod.run_hpo_annotation("s1", cfg, checkpoint)

    assert not checkpoint.done("s1", "hpo")


def test_run_hpo_annotation_checkpoints_on_real_success(tmp_path: Path, monkeypatch):
    from exomeflow.config import Config
    from exomeflow.utils import Checkpoint

    mapping = tmp_path / "genes_to_phenotype.txt"
    mapping.write_text(
        "ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\n"
        "1\tBRCA1\tHP:0003002\tBreast carcinoma\n"
    )
    monkeypatch.setattr(mod, "HPO_MAPPING_FILE", mapping)

    cfg = Config(
        input_dir=tmp_path / "fastq",
        output_dir=tmp_path / "results",
        reference=tmp_path / "ref.fa",
        dbsnp=tmp_path / "dbsnp.vcf.gz",
        mills=tmp_path / "mills.vcf.gz",
        known_indels=tmp_path / "known_indels.vcf.gz",
        annovar_bin=tmp_path / "annovar",
        annovar_db=tmp_path / "annovar" / "humandb",
    )
    cfg.vcf_dir.mkdir(parents=True, exist_ok=True)
    multianno = cfg.vcf_dir / f"s1.annovar.{cfg.annovar_buildver}_multianno.txt"
    multianno.write_text("Chr\tStart\tGene.refGene\n1\t100\tBRCA1\n")

    checkpoint = Checkpoint(cfg.checkpoint_dir)
    mod.run_hpo_annotation("s1", cfg, checkpoint)

    assert checkpoint.done("s1", "hpo")


def test_enrich_splits_multi_gene_entries(tmp_path: Path, monkeypatch):
    """
    Regression test: found via audit — ANNOVAR's Gene.refGene can be a
    composite string (';'-joined for multi-gene overlap, ','-joined with
    "(distance)" suffixes for intergenic calls), which a plain merge on the
    raw column never matches against hpo_map's single-symbol keys. Those
    rows used to silently get no HPO terms with no signal anything was missed.
    """
    mapping = tmp_path / "genes_to_phenotype.txt"
    mapping.write_text(
        "ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\n"
        "1\tBRCA1\tHP:0003002\tBreast carcinoma\n"
        "2\tTP53\tHP:0002664\tNeoplasm\n"
    )
    monkeypatch.setattr(mod, "HPO_MAPPING_FILE", mapping)

    multianno = tmp_path / "s1.annovar.hg38_multianno.txt"
    multianno.write_text(
        "Chr\tStart\tGene.refGene\n"
        "1\t100\tBRCA1;TP53\n"
        "2\t200\tBRCA1(1200),TP53(300)\n"
    )
    output = tmp_path / "s1.annovar.hpo.txt"

    assert mod.enrich("s1", multianno, output) is True

    lines = output.read_text().splitlines()
    header = lines[0].split("\t")
    hpo_idx = header.index("HPO_ID")
    for row in lines[1:]:
        cols = row.split("\t")
        assert "HP:0003002" in cols[hpo_idx] and "HP:0002664" in cols[hpo_idx]
