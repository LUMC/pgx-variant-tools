import pytest
from .context import alignment

edit_dist_data = [
    # query_seq,     ref_seq,            expected,     d_fwd, d_rev
    ("TTTAAAGCTCA", "GCCTTTAAAGCTCATC", "TTTAAAGCTCA", 5,     10), # same direction
    ("TGAGCTTTAAA", "GCCTTTAAAGCTCATC", "TTTAAAGCTCA", 10,    5),  # opposite direction
    ("TTAAGCTCA",   "GCCTTTAAAGCTCATC", "TTAAGCTCA",   7,     10), # same direction, deletions in query
    ("TGACTTTAA",   "GCCTTTAAAGCTCATC", "TTAAAGTCA",   10,    7),  # opposite direction, deletions in query
    ("TTTAAAGCTCA", "GCCTTAAGCTATC",    "TTTAAAGCTCA", 7,     9),  # same direction, deletions in ref
    ("TGAGCTTTAAA", "GCCTTAAGCTATC",    "TTTAAAGCTCA", 9,     7),  # opposite direction, deletions in ref
]

test_edit_data = [(t[0], t[1], t[3]) for t in edit_dist_data]

@pytest.mark.parametrize("query_seq,ref_seq,dist", test_edit_data)
def test_edit_distance(query_seq, ref_seq, dist):
    assert alignment.get_edit_dist(query_seq, ref_seq) == dist
    

test_strand_dist_data = [(t[0], t[1], t[3], t[4]) for t in edit_dist_data]

@pytest.mark.parametrize("query_seq,ref_seq,d_fwd, d_rev", test_strand_dist_data)
def test_strand_distance(query_seq, ref_seq, d_fwd, d_rev):
    assert alignment.get_strand_dists(query_seq, ref_seq) == (d_fwd, d_rev)


test_strand_data = [(t[0], t[1], t[2]) for t in edit_dist_data]

@pytest.mark.parametrize("query_seq,ref_seq,expected", test_strand_data)
def test_adjust_strand_dir(query_seq, ref_seq, expected):
    result = alignment.adjust_strand_dir(query_seq, ref_seq)
    assert result == expected


gap_position_data = [
    (
        [
            ("AGAGGGCCCCTAA", "AGAGGG-CCCTAA", 1),
            ("AGAGGGCCCCTAA", "AGAGGGC-CCTAA", 1),
            ("AGAGGGCCCCTAA", "AGAGGGCC-CTAA", 1),
            ("AGAGGGCCCCTAA", "AGAGGGCCC-TAA", 1) # desired
        ], 
        3
    ),
    (
        [
            ("AGAGGG-CCCTAA", "AGAGGGCCCCTAA", 1),
            ("AGAGGGC-CCTAA", "AGAGGGCCCCTAA", 1),
            ("AGAGGGCC-CTAA", "AGAGGGCCCCTAA", 1),
            ("AGAGGGCCC-TAA", "AGAGGGCCCCTAA", 1) # desired
        ], 
        3
    ),
    (
        [
            ("AGAGGGCCCCTAA", "AGAGG-CCCCTAA", 1), # desired
            ("AGAGGGCCCCTAA", "AGAG-GCCCCTAA", 1),
            ("AGAGGGCCCCTAA", "AGA-GGCCCCTAA", 1)
        ], 
        0
    ),
    (
        [
            ("AGAGGGCCCCTAA", "AGAGG--CCCTAA", 1),
            ("AGAGGGCCCCTAA", "AGAGG-C-CCTAA", 1),
            ("AGAGGGCCCCTAA", "AGAGG-CC-CTAA", 1),
            ("AGAGGGCCCCTAA", "AGAGG-CCC-TAA", 1), # desired
            ("AGAGGGCCCCTAA", "AGAG-G-CCCTAA", 1),
            ("AGAGGGCCCCTAA", "AGAG-GC-CCTAA", 1),
            ("AGAGGGCCCCTAA", "AGAG-GCC-CTAA", 1),
            ("AGAGGGCCCCTAA", "AGAG-GCCC-TAA", 1),
            ("AGAGGGCCCCTAA", "AGA-GG-CCCTAA", 1),
            ("AGAGGGCCCCTAA", "AGA-GGC-CCTAA", 1),
            ("AGAGGGCCCCTAA", "AGA-GGCC-CTAA", 1),
            ("AGAGGGCCCCTAA", "AGA-GGCCC-TAA", 1)
        ], 
        3
    ),
    (
        [
            ("TTAGAAGAAGAAGATC", "TT---AGAAGAAGATC", 1),
            ("TTAGAAGAAGAAGATC", "TTAGA---AGAAGATC", 1),
            ("TTAGAAGAAGAAGATC", "TTAGAAGA---AGATC", 1),
            ("TTAGAAGAAGAAGATC", "TTAGAAGAAGA---TC", 1) # desired
        ],
        3
    ),
    (
        [
            ("TT---AGAAGAAGATC", "TTAGAAGAAGAAGATC", 1),
            ("TTAGA---AGAAGATC", "TTAGAAGAAGAAGATC", 1),
            ("TTAGAAGA---AGATC", "TTAGAAGAAGAAGATC", 1),
            ("TTAGAAGAAGA---TC", "TTAGAAGAAGAAGATC", 1) # desired
        ],
        3
    )
]

@pytest.mark.parametrize("alignments,expected", gap_position_data)
def test_pick_leftmost(alignments, expected):
    assert alignment.pick_leftmost(alignments) == alignments[expected]


test_alignment_data = [
    ("AGAGGGCCCCTAAT", "AGAGGGCCTAAT",   "AGAGGGCCCCTAAT",  "AGAGGGCC--TAAT"),
    ("AGAGGGCTAAT", "AGAGGGCCCCTAAT",  "AGAGGGC---TAAT", "AGAGGGCCCCTAAT"),
    ("GTTAGGACTAAAAATCCG", "GTTAGGACTAAAATCCG", "GTTAGGACTAAAAATCCG", "GTTAGGACTAAAA-TCCG"),
    ("GTTAGGACTAAAAATCCG", "GTTAGGACTAAATCCG",  "GTTAGGACTAAAAATCCG", "GTTAGGACTAAA--TCCG"),
    ("GTTAGGACTAAAAATCCG", "GTTAGGACTAATCCG",   "GTTAGGACTAAAAATCCG", "GTTAGGACTAA---TCCG"),
    ("GTTAGGACTAAAAATCCG", "GTTAGGACTATCCG",    "GTTAGGACTAAAAATCCG", "GTTAGGACTA----TCCG"),
    ("GTTAGGAGGTCC", "GTTAGGTCC", "GTTAGGAGGTCC", "GTTAGG---TCC"),
    ("GTTAGGAGGAGGTCC", "GTTAGGAGGTCC", "GTTAGGAGGAGGTCC", "GTTAGGAGG---TCC"),
    ("GTTAGGAGG", "GTTAGGAGG", "GTTAGGAGG", "GTTAGGAGG"),
    ("ATGGGCTCACCAGGAAAGCAAAGACACCATGGTGGCT", "ATGGGGTCACCAGAAAGCTGACGACACGAGAGTGGCT", # representative for the exon9 conversion in *83
     "ATGGGCTCACCAGGAAAGCAAAGACACCATGGTGGCT", "ATGGGGTCACCAGAAAGCTGACGACACGAGAGTGGCT"),
    ("CCCAGCCTAATTCTTTTTTTTTTTTTTTTTTTTTTGGAGACGGAGATT", "CCCAGCCTAATTCTTTTTTTTTTTTTTTTTTTTGGAGACGGAGATT",    # varying homopolymer run with and without nearby T>C conversion
     "CCCAGCCTAATTCTTTTTTTTTTTTTTTTTTTTTTGGAGACGGAGATT", "CCCAGCCTAATTCTTTTTTTTTTTTTTTTTTTT--GGAGACGGAGATT"),
    ("CCCAGCCTAATTCTTTTTTTTTTTTTTTTTTTTTTGGAGACGGAGATT", "CCCAGCCTAATCCTTTTTTTTTTTTTTTTTTTTGGAGACGGAGATT", 
     "CCCAGCCTAATTCTTTTTTTTTTTTTTTTTTTTTTGGAGACGGAGATT", "CCCAGCCTAATCCTTTTTTTTTTTTTTTTTTTT--GGAGACGGAGATT"),
    ("CCCAGCCTAATTCTTTTTTTTTTTTTTTTTTTTGGAGACGGAGATT", "CCCAGCCTAATTCTTTTTTTTTTTTTTTTTTTTTTGGAGACGGAGATT", 
     "CCCAGCCTAATTCTTTTTTTTTTTTTTTTTTTT--GGAGACGGAGATT", "CCCAGCCTAATTCTTTTTTTTTTTTTTTTTTTTTTGGAGACGGAGATT"),
    ("CCCAGCCTAATCCTTTTTTTTTTTTTTTTTTTTGGAGACGGAGATT", "CCCAGCCTAATTCTTTTTTTTTTTTTTTTTTTTTTGGAGACGGAGATT", 
     "CCCAGCCTAATCCTTTTTTTTTTTTTTTTTTTT--GGAGACGGAGATT", "CCCAGCCTAATTCTTTTTTTTTTTTTTTTTTTTTTGGAGACGGAGATT")
]

@pytest.mark.parametrize("data", test_alignment_data)
def test_get_alignment(data):
    aln = alignment.get_alignment(data[0], data[1])
    assert aln[0] == data[2]
    assert aln[1] == data[3]


@pytest.mark.parametrize("data", test_alignment_data)
def test_alignment(data):
    aln = alignment.align(data[0], data[1], window=50)
    assert aln[0] == data[2]
    assert aln[1] == data[3]


test_find_indel_data = [
    ("AGAGGGC---TAAT", "AGAGGGCCCCTAAT", [7,8,9]),
    ("CCCAGCC-AATTCTTTTTTTTTTTTTTTTTTTTTTGGAGACGGAGAT-", "-CCAGCCTAATTCTT-TTTTTTTTTTTTTTTTTTTG-AGACGGAGATT", [0,7,15,36,47])
]

@pytest.mark.parametrize("data", test_find_indel_data)
def test_find_indels(data):
    assert list(alignment.find_indels(data[0], data[1])) == data[2]


def test_parse_cigar():
    assert alignment.parse_cigar("4=1I5X2D") == [(4,"="), (1, "I"), (5, "X"), (2, "D")]


def test_apply_cigar():
    data=("AGAGGGCTAAT", "AAGGGCCCCTTAT", "AGAGGGC---TAAT", "A-AGGGCCCCTTAT")
    cigar = "1=1I5=3D1=1X2="
    assert alignment.apply_cigar(data[0], data[1], cigar) == {"query": data[2], "target": data[3]}
