import pytest
from .context import str_search

def test_str_search_1():
    seq = "TGGCCTATATATCCCTTATAT"
    result = list(str_search.str_search("TA", seq))
    assert len(result) == 2

    assert result[0]["position"] == 5
    assert result[0]["canonical_unit"] == "TA"
    assert result[0]["num_units"] == 3
    assert result[0]["repeat_sequence"] == "TATATA"

    assert result[1]["position"] == 16
    assert result[1]["canonical_unit"] == "TA"
    assert result[1]["num_units"] == 2
    assert result[1]["repeat_sequence"] == "TATA"


def test_str_search_2():
    seq = "TGGCCTATATATCCCTTATAT"
    result = list(str_search.str_search("GCG", seq))
    assert len(result) == 0


def test_str_search_3():
    seq = "TGGCCTATATATATCCCTTATAT"
    result = list(str_search.str_search("TATA", seq))

    assert len(result) == 1
    assert result[0]["position"] == 5
    assert result[0]["canonical_unit"] == "TATA"
    assert result[0]["num_units"] == 2
    assert result[0]["repeat_sequence"] == "TATATATA"


def test_parse_hgvs_with_ref():
    hgvs = "100_101TA[8]"
    result = str_search.parse_hgvs(hgvs)

    assert result["start"] == 100
    assert result["end"] == 101
    assert result["unit"] == "TA"
    assert result["count"] == 8


def test_parsge_hgvs_no_ref():
    hgvs = "100_101[8]"
    result = str_search.parse_hgvs(hgvs)

    assert result["start"] == 100
    assert result["end"] == 101
    assert result["unit"] == ""
    assert result["count"] == 8


def test_parse_hgvs_invalid_1():
    hgvs = "100_101_TA[]"

    with pytest.raises(AssertionError):
        result = str_search.parse_hgvs(hgvs)


def test_parse_hgvs_invalid_2():
    hgvs = "100_103_TA[4]"

    with pytest.raises(AssertionError):
        result = str_search.parse_hgvs(hgvs)


def test_parse_hgvs_nonsense():
    hgvs = "gobbledygook"

    with pytest.raises(AssertionError):
        result = str_search.parse_hgvs(hgvs)