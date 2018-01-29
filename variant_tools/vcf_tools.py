"""
Utilities for vcf records
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .mutalyzer import chrom_to_accession, apply_variants
from .utils import HGVS


def generate_accession_dict(reader, build='hg19'):
    """
    Generate dictionary of chrom -> accession
    for all chromosomes in reader
    :param reader: cyvcf2 reader
    :param build: genome build
    :return: dict of chrom -> accession
    """
    return {c: chrom_to_accession(c, build) for c in reader.seqnames}


def vcf_to_hgvs(record, sample_idx, accession_dict=None):
    """
    Create delins hgvs representation of a vcf record. 
    Will select the allele(s) present in the given sample
    
    :param record: cyvcf2 record
    :param sample_idx: sample index
    :param accession_dict: optional dictionary of chrom -> accession nr

    :return: A list of hgvs objects.
    List may be empty if sample is homozygous reference for the site.
    """
    gts = record.gt_bases[sample_idx]
    if "|" in gts:
        bases = gts.split("|")
    else:
        bases = gts.split("/")
    true_alt = set(bases) - set([str(record.REF)])
    if len(str(record.REF)) == 1:
        pos_str = "{0}".format(record.POS)
    else:
        pos_str = "{0}_{1}".format(record.POS, record.POS + len(str(record.REF)) - 1)
    hgvs = []
    for alt in true_alt:
        if alt == "." or alt == "*" or alt is None:
            continue
        if accession_dict is not None:
            chrom = accession_dict[record.CHROM]
        else:
            chrom = record.CHROM
        desc_str = "{0}del{1}ins{2}".format(pos_str, record.REF, alt)
        hg = HGVS(chrom, "g", desc_str)
        hgvs.append(hg)
    return hgvs


def region_to_hgvs(reader, region, sample_idx, accession_dict=None):
    """
    Generator of hgvs descriptions for a region
    :param reader: cyvcf2 reader
    :param region: Region object
    :param sample_idx: sample index
    :param accession_dict: optional dictionary of chrom -> accession nr
    :return: generator of hgvs descriptions
    """
    region_str = "{0}:{1}-{2}".format(
        region.chromosome,
        region.start,
        region.end
    )
    for record in reader(region_str):
        for hg in vcf_to_hgvs(record, sample_idx, accession_dict):
            yield hg


def region_to_sequence(reader, region, sample_idx, accession_dict=None):
    """
    Create applied sequence from vcf region
    :param reader: cyvcf2 reader
    :param region: Region object
    :param sample_idx: sample index
    :param accession_dict: optional dictionary of chrom -> accession nr
    :return: string containing sequence
    """
    hgs = [x.description for x in region_to_hgvs(reader, region, sample_idx, accession_dict)]
    if accession_dict is not None:
        if region.chromosome not in accession_dict:
            accession_dict.update({region.chromosome: chrom_to_accession(region.chromosome)})
        return apply_variants(hgs, accession_dict[region.chromosome], region.start, region.end)
    return apply_variants(hgs, region.chromosome, region.start, region.end)


def sequence_to_fasta(sequence, name):
    """
    Create fasta string of a sequence and name
    :param sequence: string with sequence
    :param name: name of sequence
    :return: fasta string. This _will_ contain newlines
    """
    s = Seq(sequence)
    r = SeqRecord(s, id=name, description="")
    return r.format("fasta")


def region_to_fasta(reader, region, sample_idx, accession_dict=None):
    """
    Create fasta string from vcf region

    The name of the fasta sequence will be the region's value attribute.
    If this attribute is None, the string representation of the region
    will be the fasta record name.
    :param reader: cyvcf2.reader
    :param region: Region object
    :param sample_idx: sample index
    :param accession_dict: optional dictionary of chrom -> accession nr
    :return: fasta string. Contains newlines.
    """
    seq = region_to_sequence(reader, region, sample_idx, accession_dict)
    if region.value is not None:
        return sequence_to_fasta(seq, region.value)
    else:
        return sequence_to_fasta(seq, str(region))
