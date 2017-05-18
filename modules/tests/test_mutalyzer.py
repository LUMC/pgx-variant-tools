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
    (["1C>T"], "TGTTTTGGGGTGTCATATATG"),
    (["1C>T", "2G>A"], "TATTTTGGGGTGTCATATATG"),
    (["3_6del", "7_8insCC"], "CGGCCGGGTGTCATATATG")
]

@pytest.mark.parametrize("vars,expected", variant_data1)
def test_apply_variants_mutalyzer(vars, expected):
    id = mutalyzer.get_mutalyzer_region_id("NC_000022.11", 42126100, 42126120)
    assert mutalyzer.apply_variants_mutalyzer(vars, id) == expected


variant_data2 = [
    (["42126100C>T"], "TGTTTTGGGGTGTCATATATG"),
    (["42126100C>T", "42126101G>A"], "TATTTTGGGGTGTCATATATG"),
    (["42126102_42126105del"], "CGGGGGTGTCATATATG"),
    (["42126102_42126105del", "42126106_42126107insCC"], "CGGCCGGGTGTCATATATG")
]

@pytest.mark.parametrize("vars,expected", variant_data2)
def test_apply_variants(vars, expected):
    assert mutalyzer.apply_variants(vars, "NC_000022.11", 42126100, 42126120) == expected