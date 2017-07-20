"""
Generic utils
"""

from collections import namedtuple


class Region(namedtuple('Region', ['chromosome', 'start', 'end',  'value'])):
    __slots__ = ()

    def __new__(cls, chromosome, start, end, value=None):
        return super().__new__(cls, chromosome, start, end, value)

    def __str__(self):
        return "{0}:{1}-{2}".format(self.chromosome, self.start, self.end)


class HGVS(namedtuple("HGVS", ['accesssion', 'type', 'description'])):

    def __str__(self):
        return "{0}:{1}.{2}".format(self.accesssion, self.type, self.description)


def regions_from_file(inputf):
    """Generator of Region's from BED file"""
    with open(inputf) as handle:
        for line in handle:
            contents = line.strip().split("\t")
            chrom, start, end = contents[:3]
            if len(contents) == 3:
                value = None
            else:
                value = contents[3]
            yield Region(chrom, int(start), int(end), value)
