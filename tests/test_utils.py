from os.path import dirname, realpath, join

import pytest

from .context import utils


@pytest.fixture()
def bed_file():
    return join(join(dirname(realpath(__file__)), "data"), "test.bed")

hgvs_data = [
    (utils.HGVS("NC_000023.10", "g", "33038255C>A"), "NC_000023.10:g.33038255C>A"),
    (utils.HGVS("NM_004006.1", "c", "93+1G>T"), "NM_004006.1:c.93+1G>T"),
    (utils.HGVS("NC_000023.10", "g", "6775delinsGA"), "NC_000023.10:g.6775delinsGA")
]


@pytest.mark.parametrize("obj, representation", hgvs_data)
def test_hgvs(obj, representation):
    assert str(obj) == representation


region_data = [
    (utils.Region('chr1', '100', '200'), "chr1:100-200"),
    (utils.Region('chr2', '200', '300'), "chr2:200-300")
]


@pytest.mark.parametrize("obj, representation", region_data)
def test_region(obj, representation):
    assert str(obj) == representation


def test_region_from_file(bed_file):
    regions = [str(x) for x in utils.regions_from_file(bed_file)]
    assert "chr1:100-200" in regions
    assert "chr2:200-300" in regions
    assert "chr3:300-400" in regions
