"""
Utilities using the Mutalyzer web service
"""
import requests
from collections import defaultdict
import re

# cache for mutalyzer sliceChromosome results
chr_id_cache = defaultdict(str)

# url to the mutalyzer api
MUTALYZER = "https://mutalyzer.nl"


def offset_hgvs(hgvs, offset):
    """
    change the positions associated with a (simple) variant in hgvs notation
    by applying the specified offset

    :param hgvs: A string containing a variant description in hgvs notation
    :param offset: An integer offset to be applied to variant coordinates
    :return: A string containing the hgvs variant description after offsetting
    """
    if '_' in hgvs:
        # coordinate pair
        m = re.match(r'([0-9]+)[_]+([0-9]+)([a-zA-Z]+.*)', hgvs, re.M | re.I)
        assert len(m.groups()) == 3, "could not decode variant"
        return str(int(m.group(1)) + offset) + "_" + str(int(m.group(2)) + offset) + m.group(3)
    else:
        # single coordinate
        m = re.match(r'([0-9]+)([a-zA-Z]+.*)', hgvs, re.M | re.I)
        assert len(m.groups()) == 2, "could not decode variant"
        return str(int(m.group(1)) + offset) + m.group(2)


def get_mutalyzer_region_id(chrAcc, start, end):
    """
    Retrieve an ID from mutalyzer sliceChromosome SOAP service which represents
    the target region of the reference genome.
    The ID can be used as a reference against which to apply other 
    mutalyzer operations

    :param chrAcc: chromosome accession number
    :param start: start position of region (1-based inclusive)
    :param end: end position of region (1-based inclusive)
    :return: a string containing the mutalyzer id
    """

    # if this region is cached, use the cached result
    cached_id = chr_id_cache[(chrAcc, start, end)]
    if cached_id != "":
        return cached_id

    # otherwise query mutalyzer
    mutalyzer_url = "{}/json/sliceChromosome".format(MUTALYZER)
    mutalyzer_params = {
        "chromAccNo": chrAcc,
        "start": str(start),
        "end": str(end),
        "orientation": "1" 
    }

    #print("Querying mutalyzer for chromosome region")
    response = requests.get(url=mutalyzer_url, params=mutalyzer_params)
    ref_id = response.text.strip('"')
    #print("Mutalyzer returned {}".format(ref_id))
    assert ref_id.startswith("UD_"), "failed to obtain reference ID from mutalyzer"

    # cache the result
    chr_id_cache[(chrAcc, start, end)] = ref_id

    return ref_id


def apply_variants_mutalyzer(variants, chr_id):
    """
    Get the sequence obtained by applying the specified 'variants'
    to the reference sequence

    The variants to be applied must be specified with coordinates relative to
    the start of the region defined by chr_id

    :param variants: a list strings containing hgvs variant definitions
    :param chr_id: the mutalyzer id for the target reference region
    :return: A string containing the 'mutated' sequence
    """
    if len(variants) == 0:
        variants = ["1dup", "1del"] # this will give reference sequence

    mutalyzer_url = "{}/json/runMutalyzer".format(MUTALYZER)
    mutalyzer_params = {
        "variant": chr_id + ":g.[{}]".format(";".join(variants))
    }
    response = requests.get(url=mutalyzer_url, params=mutalyzer_params)
    sequence = response.json()

    return sequence["mutated"]


def apply_variants(variants, chrAcc, start, end):
    """
    Generate a 'mutated' sequence by applying the given variants to the specified genomic region

    :param variants: a list of strings containing the hgvs definitions of the variants to be applied
    :param chrAcc: Accession number for the target chromosome sequence
    :param start: start position of genomic region
    :param end: end position of genomic region
    :return: A string containing the 'mutated' sequence
    """
    chr_id = get_mutalyzer_region_id(chrAcc, start, end)
    return apply_variants_mutalyzer([offset_hgvs(v, -start + 1) for v in variants], chr_id)


def chrom_to_accession(chrom, build='hg19'):
    """
    Generate accession number from chromosome name
    :param chrom: name of chromosome
    :param build: The genome build
    :return: A string containing the accession number
    """
    mutalyzer_url = "{0}/json/chromAccession".format(MUTALYZER)
    mutalyzer_params = {
        'build': build,
        'name': chrom
    }
    response = requests.get(url=mutalyzer_url, params=mutalyzer_params)
    return response.json()
