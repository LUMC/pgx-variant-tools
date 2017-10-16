"""
Alignment of sequences
"""

from Bio import SeqIO, pairwise2
from Bio.Seq import reverse_complement
import edlib
from interval import interval
import numpy as np
import sys

# increase the max number of alignments that pairwise2 will return (default is 1000, which is insufficient when alignment is complex)
pairwise2.MAX_ALIGNMENTS = 10000


def get_edit_dist(seq1, seq2):
    """
    Get the edit distance between two sequences
    The edit distance is computed using edlib as it is considerably faster than python levenshtein
    
    :param seq1, seq2: sequences to be compared
    :type seq1, seq2:  string
    :return:           the edit distance (int)
    """
    return edlib.align(seq1, seq2, task="distance", mode="NW")["editDistance"]


def get_strand_dists(query_seq, ref_seq):
    """
    Measure the distances between a query sequence and reference sequence,
    and the reverse complement of query sequence and the reference sequence
    
    :param query_seq, ref_seq: query and reference sequences
    :type query_seq, ref_seq:  string
    :return:                   a tuple containing the pair of distances (d_fwd, d_rev) 
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
        
    :param ref_seq:   The reference sequence
    :type ref_seq:    string
    :param query_seq: The query sequence
    :type query_seq:  string
    :param verbose:   A flag indicating if results should be printed to console
    :return:          A SeqRecord containing the query sequence adjusted to match 
                      the reference direction
    """
    d_fwd, d_rev = get_strand_dists(query_seq, ref_seq)
   
    if d_rev < d_fwd:
        return reverse_complement(query_seq)

    return query_seq


def find_indels(ref_aln, qry_aln):
    """
    Identify the positions of indels in the alignment of 2 sequences
    
    :param ref_aln: An aligned sequence representing the reference
    :param qry_aln: An aligned sequence representing the query
    :yields:        The position of each gap character in the alignment
    """
    assert len(ref_aln) == len(qry_aln), "Error - Aligned sequence lengths do not match"
    for pos, pair in enumerate(zip(ref_aln, qry_aln)):
        if '-' in pair:
            yield pos
            

def pick_leftmost(alignments):
    """
    Given a list of alignments, identify the most left-justified (i.e. gaps are rightmost)
    
    To identify the rightmost positioning of gaps, the mean of the gap end positions
    for each alignment is used as a heuristic.
    
    :param alignments: A list of alignment tuples with structure (ref, query, score)
    :return:           The highest scoring left-justified alignment
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


def parse_cigar(cigar):
    """
    Parse a cigar string into a list of operations
    
    :param cigar: a cigar string
    :returns:     A list of tuples of the form (length, operation)
    """
    cigar = cigar
    ops = []
    while len(cigar) > 0:
        for i, c in enumerate(cigar):
            if c in ("=", "X", "D", "I"):
                op = cigar[:i+1]
                ops.append((int(op[:-1]), op[-1]))
                cigar = cigar[i+1:]
                break
    return ops


def apply_cigar(query, target, cigar):
    """
    Given two unaligned sequences and a cigar string, return 
    aligned sequences
    
    :param query:  query sequence
    :param target: regerence sequence
    :param cigar:  cigar string
    :returns:      a dictionary containing query and target sequences
                   with gap characters
    """
    ops = parse_cigar(cigar)
    ori1, ori2 = (query, target)
    seq1, seq2 = ("", "")
    
    for count, op in ops:
        if op in ("=", "X"):
            seq1 += ori1[:count]
            ori1 = ori1[count:]
            seq2 += ori2[:count]
            ori2 = ori2[count:]
        elif op == "I":
            seq1 += ori1[:count]
            ori1 = ori1[count:]
            seq2 += "-" * count
        elif op == "D":
            seq1 += "-" * count
            seq2 += ori2[:count]
            ori2 = ori2[count:]
    
    return {"query": seq1, "target": seq2}


default_alignment_opts = dict(
    match=3,
    mismatch=-3,
    gap_open=-9,
    gap_extend=-2,
    penalize_opening=True,
    semi_global=True
)

def get_alignment(ref_seq, query_seq, aln_opts={}):
    """
    Perform a (semi) global alignment of two sequences using bio.pairwise2
    
    :param ref_seq:   reference sequence (str)
    :param query_seq: query sequence (str)
    :param aln_opts:  additional arguments to be passed to pairwise2.align.globalms
        
    :return: The top scoring, most left-aligned alignment
    """
    opts = default_alignment_opts.copy()
    opts.update(aln_opts)
    
    alns = pairwise2.align.globalms(
        ref_seq, query_seq,
        opts["match"],
        opts["mismatch"],
        opts["gap_open"],
        opts["gap_extend"],
        penalize_end_gaps=not opts["semi_global"],
        penalize_extend_when_opening=opts["penalize_opening"]
    )
    
    # print out a warning if max_alignments is reached, as in this case correct justification is not guaranteed
    if len(alns) == pairwise2.MAX_ALIGNMENTS:
        print("Warning: MAX_ALIGNMENTS ({}) exceeded".format(pairwise2.MAX_ALIGNMENTS), file=sys.stderr)
    
    return pick_leftmost(alns)


def align(ref_seq, query_seq, window=100, aln_opts={}):
    """
    Align two sequences then perform realignment around indels
    in order to obtain consistent gap positions.
    
    1. Coarse alignment
    2. Identify the location of indels
    3. Create intervals of 'window' size around the indel locations
    4. Merge overlapping intervals
    5. Fine realignment of the sequences contained within each interval
    6. Splice the re-aligned sections back into the original alignment
    
    :param ref_seq:   The reference sequence
    :param query_seq: The query sequence
    :param window:    The number of nucleotides each side of an indel to use for re-alignment
    :param aln_opts:  Options to be passed to the alignment function
                      See variant_tools.alignment.get_alignment
    
    :returns:         A tuple containing the aligned reference sequence, the aligned query sequence
                      and the intervals used for re-alignment
    """
    
    # adjust strandedness
    query_seq = adjust_strand_dir(query_seq, ref_seq)
    
    # coarse align
    aln = edlib.align(query_seq, ref_seq, mode="NW", task="path")
    result = apply_cigar(query_seq, ref_seq, aln["cigar"])
    
    # use aligned sequences
    query_seq = result["query"]
    ref_seq = result["target"]
    
    aln_len = len(ref_seq)
    intervals = interval()
    
    # create intervals which cover the indel sites
    for pos in find_indels(ref_seq, query_seq):
        intervals |= interval[max(0, pos - window), min(aln_len, pos + window)]

    # strings to hold the patched alignments
    new_ref_seq = ""
    new_query_seq = ""
    
    # the last position covered
    last_pos = 0
    
    # re-align each interval
    for i in intervals.components:
        # start and end of the interval
        start, end = int(i[0][0]), int(i[0][1])
        
        # add the original sequences from last used position up to current start
        new_ref_seq += ref_seq[last_pos:start]
        new_query_seq += query_seq[last_pos:start]
        
        # Extract the aligned sequences within the current interval 
        # and remove any gap characters
        r = ref_seq[start:end].replace('-', '')
        q = query_seq[start:end].replace('-', '')
        
        # re-align
        aln = get_alignment(r, q, aln_opts)
        
        # splice the realigned section back into the original alignment
        new_ref_seq += aln[0]
        new_query_seq += aln[1]
        
        last_pos = end
    
    # add any remaining sequence
    new_ref_seq += ref_seq[last_pos:]
    new_query_seq += query_seq[last_pos:]
   
    return (new_ref_seq, new_query_seq, intervals)
    

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
