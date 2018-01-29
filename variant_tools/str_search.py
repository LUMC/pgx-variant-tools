import re

def parse_hgvs(hgvs):
    """
    Parse an hgvs (or hgvs-like) description of a tandem repeat unit

    :param hgvs: a string containing the variant description
    :return: a dictionary containing of the form
        {
            "start": start_position,
            "end": end_position,
            "unit": the repeat unit string
            "count": the number of repeats
        }
    """
    regex_string = "(?P<start>[0-9]+)_(?P<end>[0-9]+)(?P<unit>[A-Za-z]*)\[(?P<count>[0-9]+)\]"
    regex = re.compile(regex_string)

    match = regex.match(hgvs)

    assert match is not None, "description could not be parsed"
    if len(match.group("unit")) > 0:
        assert len(match.group("unit")) == int(match.group("end")) + 1 - int(match.group("start")), "range and sequence do not match"

    return {
        "start": int(match.group("start")),
        "end": int(match.group("end")),
        "unit": match.group("unit"),
        "count": int(match.group("count"))
    }


def str_search(repeat_unit, sequence, start_pos=0, offset=0):
    """
    Search a sequence for short tandem repeats of a canonical repeat unit.

    :param repeat_unit: The canonical repeat unit
    :param sequence:    The sequence to search within
    :param start_pos:   Position in 'sequence' at which to begin the search
    :param offset:      Optional integer offset to be added to the start position
    """
    regex_string = "({}){{2,}}".format(repeat_unit)
    regex = re.compile(regex_string)

    for result in regex.finditer(sequence, pos=start_pos):
        yield dict(
            position=result.span()[0] + offset,
            canonical_unit=repeat_unit,
            num_units=(result.span()[1] - result.span()[0]) // len(repeat_unit),
            repeat_sequence=result.group(0)
        )


def str_search_aligned(repeat_unit, sequence):
    """
    Search a (gapped) aligned sequence for short tandem repeats of a canonical repeat_unit

    :param repeat_unit: The canonical repeat unit
    :param sequence:    The (aligned) sequence to search within
    :return: a summary dict for each found repeat
    """
    regex_string = "(" + "".join(["{}[-]*".format(n) for n in repeat_unit]) + "){2,}"
    regex = re.compile(regex_string)
    for result in regex.finditer(sequence):
        gapped_seq = result.group(0)
        ungapped_seq = gapped_seq.replace("-", "")
        yield dict(
            start = result.span()[0],
            end = result.span()[1] - 1,
            canonical_unit = repeat_unit,
            num_units = len(ungapped_seq) // len(repeat_unit),
            gapped_sequence = gapped_seq,
            ungapped_sequence = ungapped_seq
        )