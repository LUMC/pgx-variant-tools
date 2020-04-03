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
        "chr10": "NC_000010.10",
        "chr11": "NC_000011.9",
        "chr12": "NC_000012.11",
        "chr13": "NC_000013.10",
        "chr14": "NC_000014.8",
        "chr15": "NC_000015.9",
        "chr16": "NC_000016.9",
        "chr17": "NC_000017.10",
        "chr18": "NC_000018.9",
        "chr19": "NC_000019.9",
        "chr2": "NC_000002.11",
        "chr20": "NC_000020.10",
        "chr21": "NC_000021.8",
        "chr22": "NC_000022.10",
        "chr3": "NC_000003.11",
        "chr4": "NC_000004.11",
        "chr5": "NC_000005.9",
        "chr6": "NC_000006.11",
        "chr8": "NC_000008.10",
        "chr7": "NC_000007.13",
        "chr9": "NC_000009.11",
        "chrM": "NC_012920.1",
        "chrX": "NC_000023.10",
        "chrY": "NC_000024.9",
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
