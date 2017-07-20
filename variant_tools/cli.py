import re

import click
import cyvcf2

from .utils import Region, regions_from_file
from .vcf_tools import generate_accession_dict, region_to_fasta


def validate_region_str(ctx, param, value):
    """
    Validate region string as a click parameter. 
    Must be of format 'chromosome:start-end'
    Returns region object
    May raise BadParameter.
    """
    if value is None:
        return None
    region_re = re.compile('^([\w\d]+):(\d+)-(\d+)$')
    match = region_re.match(value)
    if match is not None:
        return Region(match.group(1),
                      int(match.group(2)),
                      int(match.group(3)))
    else:
        raise click.BadParameter('{0} is not a '
                                 'valid region string'.format(value))


@click.command(short_help="Convert vcf records to fasta sequence")
@click.option("--vcf", "-v", type=click.Path(exists=True), required=True,
              help="Path to input VCF")
@click.option("--sample", "-s", type=click.STRING,
              help="Sample to consider. "
                   "If not given, take fist sample in VCF")
@click.option("--region", "-R", callback=validate_region_str,
              help="Region string of format chr:start-end")
@click.option('--region-file', "-L", type=click.Path(exists=True),
              help="Path to BED file containing regions")
def vcf2sequence(vcf, sample=None, region=None, region_file=None):
    """
    Convert vcf records into fasta sequences.

    Will emit a fasta record for every region in `-L` or `-R`.
    Names of sequences will be taken from either:
        - an optional fourth column in the BED file
        - the string representation of the region
    """
    if region is None and region_file is None:
        raise ValueError("Either -R or -L must be set")

    reader = cyvcf2.VCF(vcf)
    if sample is not None:
        sample_idx = reader.samples.index(sample)
    else:
        sample_idx = 0

    acdict = generate_accession_dict(reader)

    if region:
        print(region_to_fasta(reader, region, sample_idx, acdict), end='')

    if region_file:
        for reg in regions_from_file(region_file):
            print(region_to_fasta(reader, reg, sample_idx, acdict), end='')
