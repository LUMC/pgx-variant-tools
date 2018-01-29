"""
Calling of variants using a minimal set of HGVS-like descriptions

The minimal notation uses three operators:
    substitution
    insertion
    deletion
"""
from enum import Enum


class Opcodes(Enum):
    """
    Operator codes
    """
    EQ  = 0
    SUB = 1
    INS = 2
    DEL = 3
    SKIP = 4


class Operator(object):
    """
    base class for all operators
    """
    opcode = None
    
    def __init__(self, pos=0, ref=None, alt=None):
        self.pos = pos
        
        self.ref = []
        if ref is not None:
            self.ref.append(ref)
            
        self.alt = []
        if alt is not None:
            self.alt.append(alt)

    def __str__(self):
        raise NotImplementedError()


class Substitution(Operator):
    """
    substitution operator
    """
    opcode = Opcodes.SUB
    
    def __str__(self):
        return "{}{}>{}".format(self.pos, self.ref[0], self.alt[0])
    
    
class Deletion(Operator):
    """
    deletion operator
    """
    opcode = Opcodes.DEL
    
    def __str__(self):
        if len(self.ref) == 1:
            return "{}del{}".format(self.pos, self.ref[0])
        
        return "{}_{}del{}".format(self.pos, self.pos + len(self.ref) - 1, "".join(self.ref))


class Insertion(Operator):
    """
    insertion operator
    """
    opcode = Opcodes.INS
    
    def __str__(self):
        return "{}_{}ins{}".format(self.pos, self.pos + 1, "".join(self.alt))

  
def get_opcode(ref_base, alt_base):
    """
    Get the opcode for the operation needed to convert a the given reference base into the given alt base
    
    :param ref_base: Character representing the reference base
    :param alt_base: Character representing the alternative base
    :return: 
    """
    if ref_base == alt_base:
        return Opcodes.EQ
    if ref_base == "$" or alt_base == "$":
        return Opcodes.SKIP
    if ref_base == "-":
        return Opcodes.INS
    if alt_base == "-":
        return Opcodes.DEL
    return Opcodes.SUB


def call_variants(alignment, offset=0):
    """
    Describe an alignment as a list of simple operations,
    where an operation is one of Substitution, Insertion or Deletion
    
    :param alignment: An alignment object tuple of the structure (ref_seq_aligned, alt_seq_aligned)
    :param offset:    An integer offset to be applied to all positions
    :return:          A list of variant operators. 
                      For compatability with HGVS, operators use a 1-based coordinate system
    """
    ops = []                # list of variant operations
    cur_op = None           # the current variant operation
    s_ref = alignment[0]    # ref sequence
    s_alt = alignment[1]    # alt sequence
    r_pos = -1              # current position on reference sequence
    
    # iterate over the alignment
    for i in range(len(alignment[0])):
        c_ref = s_ref[i]         # ref character at position i
        c_alt = s_alt[i]         # alt character at position i
        
        # increment ref sequence position
        if c_ref != "-":
            r_pos += 1
        
        # determine the opcode for ref vs alt at current position
        op = get_opcode(c_ref, c_alt)
        
        # store current op if new op is different
        if cur_op is not None and op != cur_op.opcode:
            ops.append(cur_op)
            cur_op = None
            
        # substitution: single-nucleotide operation
        # create Substitution operator, add it to ops
        if op == Opcodes.SUB:
            ops.append(Substitution(r_pos + 1 + offset, c_ref, c_alt))
        
        # insertion: multi-nucleotide operation,
        # create Insertion operator and hold onto it for further work
        if op == Opcodes.INS:
            if cur_op is None:
                cur_op = Insertion(r_pos + 1 + offset, alt=c_alt)
            else:
                cur_op.alt.append(c_alt)
        
        # deletion: multi-nucleotide operation,
        # create Deletion operator and hold onto it for further work
        if op == Opcodes.DEL:
            if cur_op is None:
                cur_op = Deletion(r_pos + 1 + offset, ref=c_ref)
            else:
                cur_op.ref.append(c_ref)
    
    # append the last operation
    if cur_op is not None:
        ops.append(cur_op)
    
    return ops
