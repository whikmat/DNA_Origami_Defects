from nupack import *

def test_version():
    import nupack
    assert '@' not in nupack.__version__

################################################################################

def test_strands():
    A = RawStrand('AGTCTAGGATTCGGCGTGGGTTAA')
    B = RawStrand('TTAACCCACGCCGAATCCTAGACTCAAAGTAGTCTAGGATTCGGCGTG')
    C = RawStrand('AGTCTAGGATTCGGCGTGGGTTAACACGCCGAATCCTAGACTACTTTG')

    a1 = Domain('AGTCTAGGATTCGGCGT', name='a1')
    a2 = Domain('GGGTTAA', name='a2')
    D = TargetStrand([a1, a2], name='D') # mostly useful in a design context


    c1 = RawComplex([A])
    c2 = RawComplex([A, B, B, C])
    c3 = RawComplex([A, A])
    c4 = RawComplex([A, B, C])
    c1a = RawComplex([A, B])
    c1b = RawComplex([B, A])
    c1c = RawComplex([A, B])

    assert c1a == c1b
    assert c1a == c1a
    assert c1a == c1c

################################################################################