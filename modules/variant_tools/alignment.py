"""
Alignment of sequences
"""

from Bio import SeqIO, pairwise2
from Bio.Seq import reverse_complement
import edlib
import numpy as np
import sys

# increase the max number of alignments that pairwise2 will return (default is 1000, which is insufficient when alignment is complex)
pairwise2.MAX_ALIGNMENTS = 10000


def get_edit_dist(seq1, seq2):
    """
    Get the edit distance between two sequences
    The edit distance is computed using edlib as it is considerably faster than python levenshtein
    
    :param seq1, seq2: sequences to be compared
    :type seq1, seq2: string
    :return: the edit distance (int)
    """
    return edlib.align(seq1, seq2, task="distance", mode="NW")["editDistance"]


def get_strand_dists(query_seq, ref_seq):
    """
    Measure the distances between a query sequence and reference sequence,
    and the reverse complement of query sequence and the reference sequence
    
    :param query_seq, ref_seq: query and reference sequences
    :type query_seq, ref_seq: string
    :return: a tuple containing the pair of distances (d_fwd, d_rev) 
    """
    fwd = query_seq
    rev = reverse_complement(fwd)
    d_fwd = get_edit_dist(fwd, ref_seq)
    d_rev = get_edit_dist(rev, ref_seq)
    return (d_fwd, d_rev)
    
    
def adjust_strand_dir(query_seq, ref_seq):
    """
    Adjust the strandedness of 'query_seq' to match that of 'ref_seq'
    This is achieved by determining the edit distance between query_seq and its reverse_complement
    to 'ref_seq', and choosing the direction which gives the smallest distance
        
    :param ref_seq: A biopython SeqRecord containing the reference sequence
    :param query_seq: A biopython SeqRecord containing the query sequence
    :param verbose: A flag indicating if results should be printed to console
    :return: A SeqRecord containing the query sequence adjusted to match the reference direction
    """
    d_fwd, d_rev = get_strand_dists(query_seq, ref_seq)
   
    if d_rev < d_fwd:
        return reverse_complement(query_seq)

    return query_seq


def pick_leftmost(alignments):
    """
    Given a list of alignments, identify the most left-justified (i.e. gaps are rightmost)
    
    To identify the rightmost positioning of gaps, the mean of the gap end positions
    for each alignment is used as a heuristic.
    
    :param alignments: A list of alignment tuples with structure (ref, query, score)
    :return: The highest scoring left-justified alignment
    """
    best_score = 0
    best_alignment = alignments[0]
    for aln in alignments:
        gaps = []
        gap = None
        for i in range(len(aln[0])):
            if aln[0][i] == '-' or aln[1][i] == '-':
                if gap is None:
                    gap = []
                    gap.append(i) # gap start
            else:
                if gap is not None:
                    gap.append(i-1) # gap end
                    gaps.append(gap)
                    gap = None
        
        if gap is not None:
            gap.append(i-1) # gap end
            gaps.append(gap)
        
        if len(gaps) == 0:
            break
            
        score = np.mean([g[-1] for g in gaps])
        if score > best_score:
            best_score = score
            best_alignment = aln
           
    return best_alignment


def get_alignment(ref_seq, query_seq, match=3, mismatch=-3, gap_open=-9, gap_extend=-2, semi_global=False):
    """
    Perform a (semi) global alignment of two sequences
    
    :param ref_seq: reference sequence (str)
    :param query_seq: query sequence (str)
    :param match: match score
    :param mismatch: mismatch penalty
    :param gap_open: gap opening penalty
    :param gap_extend: gap extension penalty
    :param semi_global: boolean flag indicating that a semi-global alignment (i.e. terminal gaps are not scored) should be performed
                        Default is False, resulting in a full global alignment
    
    :return: The top scoring, most left-aligned alignment
    """
    alns = pairwise2.align.globalms(
        ref_seq, query_seq,
        match, mismatch, gap_open, gap_extend,
        penalize_end_gaps=semi_global,
        penalize_extend_when_opening = True
    )
    
    if len(alns) == pairwise2.MAX_ALIGNMENTS:
        print("Warning: MAX_ALIGNMENTS ({}) exceeded".format(pairwise2.MAX_ALIGNMENTS), file=sys.stderr)
    
    return pick_leftmost(alns)


def pretty_alignment(seq1, seq2, start=0, end=None, width=200, match=" ", mismatch="*", gap="#"):
    """
    Pretty string representation of the alignment
    :param seq1, seq2: aligned sequence strings
    :param start:      the alignment position at which to start
    :param end:        the alignment position at which to end
    :param width:      the number of characters per line
    :param match:      character used to represent a match
    :param mismatch:   character used to represent a mismatch
    :param gap:        character used to represent a gap
    """
    pos = start
    if end is None:
        end = len(seq1)
    
    def match_symbol(base1, base2):
        if base1 == base2:
            return match
        if base1 == "-" or base2 == "-":
            return gap
        return mismatch
        
    pretty = []

    while pos < end:
        num_bases = width if pos + width < end else end-pos
        pretty.append(seq1[pos: pos + num_bases])
        pretty.append("".join([match_symbol(seq1[pos+i], seq2[pos+i]) for i in range(num_bases)]))
        pretty.append(seq2[pos:pos+num_bases])
        pretty.append("")
        pos += width
        
    return "\n".join(pretty)
