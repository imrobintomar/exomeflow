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

    mod.enrich("s1", multianno, output)

    lines = output.read_text().splitlines()
    header = lines[0].split("\t")
    assert "HPO_ID" in header and "HPO_terms" in header

    brca1_row = next(row for row in lines[1:] if "BRCA1" in row)
    assert "HP:0003002" in brca1_row and "HP:0100615" in brca1_row

    unknown_row = next(row for row in lines[1:] if "UNKNOWN_GENE" in row)
    # pandas writes NaN for unmatched merges
    assert unknown_row.split("\t")[-1] in ("", "nan")


def test_enrich_skips_gracefully_when_mapping_missing(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(mod, "HPO_MAPPING_FILE", tmp_path / "does_not_exist.txt")

    multianno = tmp_path / "s1.annovar.hg38_multianno.txt"
    multianno.write_text("Chr\tStart\tGene.refGene\n1\t100\tBRCA1\n")
    output = tmp_path / "s1.annovar.hpo.txt"

    mod.enrich("s1", multianno, output)  # must not raise
    assert output.exists()
    assert "HPO_ID" not in output.read_text().splitlines()[0]


def test_enrich_skips_when_multianno_missing(tmp_path: Path):
    output = tmp_path / "out.txt"
    mod.enrich("s1", tmp_path / "missing.txt", output)
    assert not output.exists()
