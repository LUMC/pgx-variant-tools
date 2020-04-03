import pytest
from .context import mutalyzer

offset_data = [
    ("999C>T", 100, "1099C>T"),
    ("1234_1236delinsTTTCT", -10, "1224_1226delinsTTTCT"),
    ("1000_1001insGGG", -500, "500_501insGGG")
]

@pytest.mark.parametrize("var,offset,expected", offset_data)
def test_offset_hgvs(var, offset, expected):
    assert mutalyzer.offset_hgvs(var, offset) == expected


def test_get_mutalyzer_region_id():
    assert mutalyzer.get_mutalyzer_region_id("NC_000022.11", 42126062, 42132599) == "UD_143620441442"


variant_data1 = [
    ([], "CGTTTTGGGGTGTCATATATG"),
    (["1C>T"], "TGTTTTGGGGTGTCATATATG"),
    (["1C>T", "2G>A"], "TATTTTGGGGTGTCATATATG"),
    (["3_6del", "7_8insCC"], "CGGCCGGGTGTCATATATG")
]

@pytest.mark.parametrize("vars,expected", variant_data1)
def test_apply_variants_mutalyzer(vars, expected):
    id = mutalyzer.get_mutalyzer_region_id("NC_000022.11", 42126100, 42126120)
    assert mutalyzer.apply_variants_mutalyzer(vars, id) == expected


variant_data2 = [
    ([], "CGTTTTGGGGTGTCATATATG"),
    (["42126100C>T"], "TGTTTTGGGGTGTCATATATG"),
    (["42126100C>T", "42126101G>A"], "TATTTTGGGGTGTCATATATG"),
    (["42126102_42126105del"], "CGGGGGTGTCATATATG"),
    (["42126102_42126105del", "42126106_42126107insCC"], "CGGCCGGGTGTCATATATG")
]

@pytest.mark.parametrize("vars,expected", variant_data2)
def test_apply_variants(vars, expected):
    assert mutalyzer.apply_variants(vars, "NC_000022.11", 42126100, 42126120) == expected


chrom_data = [
    ('chr1', "NC_000001.10"),
    ('chr2', "NC_000002.11"),
    ('chr3', "NC_000003.11"),
    ('chr4', "NC_000004.11"),
    ('chr5', "NC_000005.9"),
    ('chr6', "NC_000006.11"),
    ('chr7', "NC_000007.13"),
    ('chr8', "NC_000008.10"),
    ('chr9', "NC_000009.11"),
    ('chr10', "NC_000010.10"),
    ('chr11', "NC_000011.9"),
    ('chr12', "NC_000012.11"),
    ('chr13', "NC_000013.10"),
    ('chr14', "NC_000014.8"),
    ('chr15', "NC_000015.9"),
    ('chr16', "NC_000016.9"),
    ('chr17', "NC_000017.10"),
    ('chr18', "NC_000018.9"),
    ('chr19', "NC_000019.9"),
    ('chr20', "NC_000020.10"),
    ('chr21', "NC_000021.8"),
    ('chr22', "NC_000022.10"),
    ('chrX', "NC_000023.10"),
    ('chrY', "NC_000024.9")
]


@pytest.mark.parametrize("chrom, acc", chrom_data)
def test_chrom_to_accession(chrom, acc):
    assert mutalyzer.chrom_to_accession(chrom) == acc


no_chrom_data = [
    ('chrUn_gl000239', "NT_167233.1")
]


@pytest.mark.parametrize("chrom, acc", no_chrom_data)
def test_chrom_to_accession_fail(chrom, acc):
    try:
        ret = mutalyzer.chrom_to_accession(chrom)
        assert ret != acc
    except Exception as e:
        pass
