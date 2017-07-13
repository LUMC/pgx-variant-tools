from os.path import join, dirname, realpath

from click import BadParameter
from click.testing import CliRunner
import pytest

from variant_tools.cli import validate_region_str, vcf2sequence
from variant_tools.utils import Region


@pytest.fixture()
def vcf():
    return join(join(dirname(realpath(__file__)), "data"), "mini.vcf.gz")


@pytest.fixture()
def bed():
    return join(join(dirname(realpath(__file__)), "data"), "test2.bed")


def test_validate_region_str():
    assert validate_region_str(None, None, "chr1:100-200") == Region("chr1", 100,200)
    with pytest.raises(BadParameter):
        validate_region_str(None, None, "wontparse")


def test_single_region_vcf2sequence(vcf):
    runner = CliRunner()
    result = runner.invoke(vcf2sequence, ["-v", vcf, "-R", 'chr1:14500-14600'])
    assert result.exit_code == 0


def test_from_file_vcf2sequence(vcf, bed):
    runner = CliRunner()
    result = runner.invoke(vcf2sequence, ["-v", vcf, "-L", bed])
    assert result.exit_code == 0

