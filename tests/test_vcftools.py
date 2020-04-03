from os.path import join, dirname, realpath

import cyvcf2
import pytest

from .context import vcf_tools, utils


@pytest.fixture()
def vcf():
    return cyvcf2.VCF(join(join(dirname(realpath(__file__)), "data"), "mini.vcf.gz"))


@pytest.fixture()
def acdict():
    return {"chr1": "NC_000001.10"}


def test_accession_dict(vcf):
    assert vcf_tools.generate_accession_dict(vcf) == {
        "chr1": "NC_000001.10",
    }


def test_vcf_to_hgvs(vcf, acdict):
    first = next(vcf)
    assert str(vcf_tools.vcf_to_hgvs(first, 0, acdict)[0]) == "NC_000001.10:g.14574delAinsG"
    second = next(vcf)
    assert str(vcf_tools.vcf_to_hgvs(second, 0, acdict)[0]) == "NC_000001.10:g.14590delGinsA"
    third = next(vcf)
    assert str(vcf_tools.vcf_to_hgvs(third, 0, acdict)[0]) == "NC_000001.10:g.14599delTinsA"


def test_region_to_hgvs(vcf, acdict):
    reg = utils.Region('chr1', 14500, 14600)
    assert sorted(list(map(str, vcf_tools.region_to_hgvs(vcf, reg, 0, acdict)))) == sorted([
        "NC_000001.10:g.14574delAinsG",
        "NC_000001.10:g.14590delGinsA",
        "NC_000001.10:g.14599delTinsA"
    ])


def test_region_to_sequence(vcf, acdict):
    reg = utils.Region('chr1', 14500, 14600)
    assert vcf_tools.region_to_sequence(vcf, reg, 0, acdict) == "GCACAGGCAGACAGAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTTCCGCTCCTTGAAGCTGGTCTCCGCACAGTGCTGGTTCCATCACCCCCAC"


def test_sequence_to_fasta():
    seq = "AAAAAA"
    name = "test1"
    assert vcf_tools.sequence_to_fasta(seq, name) == ">test1\nAAAAAA\n"


def test_region_to_fasta(vcf, acdict):
    reg = utils.Region('chr1', 14500, 14600)
    assert vcf_tools.region_to_fasta(vcf, reg, 0, acdict) == ">chr1:14500-14600\nGCACAGGCAGACAGAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTTCCGCTCCT\nTGAAGCTGGTCTCCGCACAGTGCTGGTTCCATCACCCCCAC\n"
    reg2 = utils.Region('chr1', 14500, 14600, "test")
    assert vcf_tools.region_to_fasta(vcf, reg2, 0, acdict) == ">test\nGCACAGGCAGACAGAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTTCCGCTCCT\nTGAAGCTGGTCTCCGCACAGTGCTGGTTCCATCACCCCCAC\n"
