import pytest
from .context import variant_calling

def test_Operator():
    o = variant_calling.Operator()
    assert o.pos == 0 
    assert o.ref == []
    assert o.alt == [] 
    assert o.opcode is None
    try:
        s = str(o)
    except NotImplementedError:
        pass


def test_Substitution():
    s = variant_calling.Substitution()
    assert s.pos == 0 
    assert s.ref == []
    assert s.alt == []
    assert s.opcode == variant_calling.Opcodes.SUB
    
    s = variant_calling.Substitution(10, "A", "T")
    assert str(s) == "10A>T"


def test_Deletion():
    d = variant_calling.Deletion()
    assert d.pos == 0
    assert d.ref == []
    assert d.alt == []
    assert d.opcode == variant_calling.Opcodes.DEL
    
    d.pos = 10
    d.ref = ["T"]
    d.alt = ["-"]
    assert str(d) == "10delT"
    
    d.ref += ["C", "G"]
    d.alt += ["-", "-"]
    assert str(d) == "10_12delTCG"


def test_Insertion():
    i = variant_calling.Insertion()
    assert i.pos == 0 
    assert i.ref == []
    assert i.alt == []
    assert i.opcode == variant_calling.Opcodes.INS
    
    i.pos = 10
    i.ref = ["-"]
    i.alt = ["T"]
    assert str(i) == "10_11insT"
    
    i.ref += ["-", "-"]
    i.alt += ["C", "G"]
    assert str(i) == "10_11insTCG"


def test_get_opcode():
    assert variant_calling.get_opcode("A", "A") == variant_calling.Opcodes.EQ
    assert variant_calling.get_opcode("-", "A") == variant_calling.Opcodes.INS
    assert variant_calling.get_opcode("A", "-") == variant_calling.Opcodes.DEL
    assert variant_calling.get_opcode("A", "G") == variant_calling.Opcodes.SUB
    assert variant_calling.get_opcode("A", "$") == variant_calling.Opcodes.SKIP
    assert variant_calling.get_opcode("$", "A") == variant_calling.Opcodes.SKIP


variant_test_data = [
    #  ref              query            offset  result
    [("TTTCTTT",       "TTT-TTT"),       0,     ["4delC"]],
    [("TTTCTTT",       "TTTATTT"),       0,     ["4C>A"]],
    [("TTT-TTT",       "TTTATTT"),       0,     ["3_4insA"]],
    [("TTTCTTT",       "TTT-TTT"),       2500,  ["2504delC"]],
    [("TTTCTTT",       "TTTATTT"),       2500,  ["2504C>A"]],
    [("TTT-TTT",       "TTTATTT"),       2500,  ["2503_2504insA"]],
    [("AGAGGGCCCCTAA", "AGAGG-CCC-TAA"), 0,     ["6delG", "10delC"]],
    [("AGAGGGCCC-TAA", "AGAGG-CCCCTAA"), 0,     ["6delG", "9_10insC"]],
    [("AGAGGGCCCCTAA", "AGAGG-CCC-TAA"), 10000, ["10006delG", "10010delC"]],
    [("AGAGGGCCC-TAA", "AGAGG-CCCCTAA"), 10000, ["10006delG", "10009_10010insC"]],
    [("AGAGGGCCCCTAA", "AGAG----CCTAA"), 0,     ["5_8delGGCC"]],
    [("AGAGGGCCCCTAA", "AGAG----CCTAA"), 9090,  ["9095_9098delGGCC"]],
    [("AGAG----CCTAA", "AGAGGGCCCCTAA"), 0,     ["4_5insGGCC"]],
    [("AGAG----CCTAA", "AGAGGGCCCCTAA"), 770,   ["774_775insGGCC"]],
    [("AGAGGGCCCCTAA", "AGAGGTCACGTAA"), 0,     ["6G>T", "8C>A", "10C>G"]],
    [("AGAGGGCCCCTAA", "AGAGGTCACGTAA"), 9999,  ["10005G>T", "10007C>A", "10009C>G"]],
    [("AGAGGGCCCCTAA", "AGAGGTCACG$$A"), 9999,  ["10005G>T", "10007C>A", "10009C>G"]],
    [("AGAGGGCCCCTAA", "AGAG----$$TAA"), 0,     ["5_8delGGCC"]],
    [("$$AGGGCCCCTAA", "AGAG----$$TAA"), 0,     ["5_8delGGCC"]],
]


@pytest.mark.parametrize("sequences,offset,expected", variant_test_data)
def test_call_variants(sequences, offset, expected):
    variants = [str(v) for v in variant_calling.call_variants(sequences, offset)]
    assert variants == expected

