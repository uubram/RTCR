#!/usr/bin/env python
#
# Description: tests RTCR pipeline.

import pytest
import tempfile
import os
from contextlib import contextmanager

import logging
import random
# Some modules use logging
logging.basicConfig(level = logging.DEBUG)

@contextmanager
def tempdata(data):
    f = tempfile.NamedTemporaryFile(delete = False)
    f.write(data)
    f.close()
    yield f.name
    os.unlink(f.name)

def test_fasta_io():
    from rtcr.fileio import zopen, FastaFormat 

    with tempdata("\n".join([
        "Comment line of fasta",
        "Another comment line",
        ">seq1 this is the first sequence",
        "ACCCTTGGGCCGGGAAATT",
        ">seq2\tstretching multiple lines",
        "accttaaaggaacc",
        "aaccttaaaaaaaa",
        ">seq3 contains upper and lower case",
        "ACCTAAAAaaaaggaaAAACCCGGT"])) as fn:
        
        f = FastaFormat.records_in(zopen(fn, 'r'))
        s = f.next()
        assert(s.name == "seq1 this is the first sequence")
        assert(s.seq == "ACCCTTGGGCCGGGAAATT")

        s = f.next()
        assert(s.name == "seq2\tstretching multiple lines")
        assert(s.seq == "accttaaaggaaccaaccttaaaaaaaa".upper())

        s = f.next()
        assert(s.name == "seq3 contains upper and lower case")
        assert(s.seq == "ACCTAAAAaaaaggaaAAACCCGGT".upper())

    #TODO: test faulty fasta files (also for wrong alphabet?)

def test_fastq_io():
    from rtcr.fileio import zopen, FastqFormat 

    with tempdata("\n".join([
            "@seq1",
            "A"*42,
            "+",
            "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ",
            "@seq2 upper and lower case",
            "AACCTTAAaaagaaagaaaTAaaaa",
            "+",
            "JJJJJJJJJJJJJJJJJJJJJJJJJ",
            "@seq3 out of range Phred quality",
            "ACTG",
            "+",
            "JKJJ"])) as fn:
        f = FastqFormat.records_in(zopen(fn, 'r'))
        s = f.next()
        
        assert(s.name == "seq1")
        assert(s.qual == tuple(range(42)))
        
        s = f.next()
        assert(s.name == "seq2 upper and lower case")
        assert(s.seq == "AACCTTAAaaagaaagaaaTAaaaa".upper())
        
        with pytest.raises(Exception) as excinfo:
            s = f.next()
        
        assert(excinfo.value.message == "Quality score (42) not in [0, 41]")

        #TODO: test for fastq files with Solexa and Illumina1.3
        #TODO: test faulty fastq files (also for wrong alphabet?)

def test_junctions2cloneset():
    from rtcr.fileio import zopen, junctions2cloneset

    jncrec = lambda read_name, nt_jnc, q_jnc:\
            "\t".join([
                read_name,
                "None", #vid
                "None", #jid
                nt_jnc,
                "".join(map(lambda q:chr(q + 33), q_jnc)),
                "0", #jnc_ve
                "0", #jnc_js
                "", #nt_jnc_gr
                "", #q_jnc_gr
                "0", #jnc_ve_gr
                "0", #jnc_js_gr
                "0", #v_score
                "0" #j_score
                ])
    header = "read_name\tvid\tjid\tnt_jnc\tq_jnc\tjnc_ve\tjnc_js\t\
nt_jnc_gr\tq_jnc_gr\tjnc_ve_gr\tjnc_js_gr\tv_score\tj_score"
    with tempdata("\n".join([header,
            jncrec("read1", "ACTG", [40]*4),
            jncrec("read2", "AAA", [10, 20, 30]),
            jncrec("read3", "ACTG", [30, 30, 41, 30]),
            jncrec("read4", "AAA", [15, 15, 41]),
            jncrec("read5", "ACCCA", [10]*5),
            jncrec("read6", "CTGA", [30]*4)])) as fn:
        cloneset = junctions2cloneset(zopen(fn, 'r'), 33)
        assert cloneset.count == 4
        assert cloneset.sequence_count == 6
        assert cloneset.base_count == 23 
        assert set([(clone.seq, tuple(clone.qual), clone.count) \
                for clone in cloneset]) == set([
                ("ACTG", (40, 40, 41, 40), 2),
                ("AAA", (15, 20, 41), 2),
                ("ACCCA", (10,)*5, 1),
                ("CTGA", (30,)*4, 1)])

def test_QSequence():
    from rtcr.seq import QSequence

    s1 = QSequence(name = "seq1", seq = "ACCCT", qual = [40]*5)
    assert s1.name == "seq1"
    assert s1.seq == "ACCCT"
    assert s1.qual == (40,)*5

    s1a = QSequence("", seq = "ACCCT", qual = [40]*5)
    assert s1a.name != s1.name
    assert s1a.seq == s1.seq
    assert s1a.qual == s1.qual
    assert s1a != s1 and not s1a == s1

    s1b = QSequence("seq1", seq = "GGGGG", qual = [40]*5)
    assert s1b.name == s1.name
    assert s1b.seq != s1.name
    assert s1b.qual == s1.qual
    assert s1b != s1 and not s1b == s1

    s1c = QSequence("seq1", seq = "ACCCT", qual = [10]*5)
    assert s1c.name == s1.name
    assert s1c.seq == s1.seq
    assert s1c.qual != s1.qual
    assert s1c != s1 and not s1c == s1

    s1d = QSequence("seq1", seq = "ACCCT", qual = [40]*5)
    assert s1d.name == s1.name
    assert s1d.seq == s1.seq
    assert s1d.qual == s1.qual
    assert s1d == s1 and not s1d != s1

    # Test reverse complement
    s2 = QSequence("", seq = "AACTGNACT", qual = range(9))
    s2rc = s2.reverse_complement()
    assert s2.name == ""
    assert s2.seq == "AACTGNACT"
    assert s2.qual == tuple(range(9))
    assert s2rc.name == ""
    assert s2rc.seq == "AGTNCAGTT"
    assert s2rc.qual == tuple(range(9)[::-1])
    assert s2rc != s2 and not s2rc == s2
    s2rc = s2.reverse_complement(name = "abc")
    assert s2rc.name == "abc"
    assert s2rc.seq == "AGTNCAGTT"

    # Test forcing Phred quality of ambiguous bases ("N") to 0
    s3 = QSequence("", seq = "NANCN", qual = [5,40,0,40,2],
            force_N_to_zero = True)
    assert s3.name == ""
    assert s3.seq == "NANCN"
    assert s3.qual == (0,40,0,40,0)

def test_PfmElement():
    import pickle
    from rtcr.seq.seq import PfmElement

    el1 = PfmElement()
    assert el1.count == 0
    assert el1.mutation_count == 0
    assert el1.base == "N" 
    assert el1.q == 0

    el1.inc_freq(base = "A", q = 30, count = 10)
    assert el1.count == 10
    assert el1.mutation_count == 0
    assert el1.base == "A"
    assert el1.q == 30

    el1.inc_freq(base = "A", q = 30, count = 20)
    assert el1.count == 30
    assert el1.mutation_count == 0
    assert el1.base == "A"
    assert el1.q == 30

    el1.inc_freq(base = "T", q = 30, count = 25)
    assert el1.count == 55
    assert el1.mutation_count == 25 
    assert el1.base == "A"
    assert el1.q == 30

    el1.inc_freq(base = "T", q = 20, count = 5)
    assert el1.count == 60
    assert el1.mutation_count == 30 
    assert el1.base == "N"
    assert el1.q == 0

    el1.inc_freq(base = "T", q = 31, count = 1)
    assert el1.count == 61
    assert el1.mutation_count == 30
    assert el1.base == "T"
    assert el1.q == 31

    # Test adding of another PfmElement
    el2 = PfmElement()
    el2.inc_freq("G", q = 40, count = 31)
    el1.add(el2)
    assert el1.count == 92
    assert el1.mutation_count == 61 
    assert el1.base == "G" # base should win based on higher quality score
    assert el1.q == 40 
    
    el3 = PfmElement()
    el3.inc_freq("A", 30, 10)
    el3.inc_freq("C", 25, 10)
    el3.inc_freq("T", 29, 40)
    el1.add(el3)
    assert el1.count == 152
    assert el1.mutation_count == 81 
    assert el1.base == "T"
    assert el1.q == 31

    # Test adding unresolved bases
    el4 = PfmElement()
    el4.inc_freq("N", q = 0, count = 10)
    assert el4.count == 10
    assert el4.mutation_count == 0
    assert el4.base == "N"
    assert el4.q == 0

    el4.inc_freq("T", q = 30, count = 2)
    assert el4.count == 12
    assert el4.mutation_count == 10
    assert el4.base == "T"
    assert el4.q == 30

    el4.inc_freq("C", q = 30, count = 2)
    assert el4.count == 14
    assert el4.mutation_count == 12
    assert el4.base == "N"
    assert el4.q == 0

    # Explicitly test resolving ties on frequency and base quality
    el5 = PfmElement()
    el5.inc_freq("A", q = 10, count = 2)
    el5.inc_freq("T", q = 10, count = 2)
    assert el5.count == 4
    assert el5.mutation_count == 2
    assert el5.base == "N"
    assert el5.q == 0
    
    el5.inc_freq("C", q = 11, count = 2)
    assert el5.count == 6
    assert el5.mutation_count == 4
    assert el5.base == "C"
    assert el5.q == 11

    el5.inc_freq("G", q = 11, count = 2)
    assert el5.count == 8
    assert el5.mutation_count == 6
    assert el5.base == "N"
    assert el5.q == 0 

    el5.inc_freq("T", q = 5, count = 1)
    assert el5.count == 9
    assert el5.mutation_count == 6
    assert el5.base == "T"
    assert el5.q == 10

    # Test pickling
    for el in [el1, el2, el3, el4, el5]:
        pickle.dumps(el)

def helperfunc_test_ConsensusQSequence(load_c):
    import pickle
    from itertools import permutations
    from rtcr.seq import QSequence
    if load_c:
        from cseq import ConsensusQSequence
    else:
        from rtcr.seq.seq import ConsensusQSequence

    # Test various ways of initialization
    c = ConsensusQSequence()
    assert c.count == 0
    assert c.base_count == 0
    assert c.mutation_count == 0
    assert c.qual == [] 
    assert c.seq == ""

    c = ConsensusQSequence(("ACCTG", [1,2,3,4,5]))
    assert c.count == 1
    assert c.base_count == 5
    assert c.mutation_count == 0
    assert c.qual == [1,2,3,4,5]
    assert c.seq == "ACCTG"

    c = ConsensusQSequence(("GGCGGT", [1,2,3,4,5,6], 7))
    assert c.count == 7
    assert c.base_count == 6 * 7
    assert c.mutation_count == 0
    assert c.qual == [1,2,3,4,5,6]
    assert c.seq == "GGCGGT"

    c = ConsensusQSequence(QSequence("", "ACA", [1,1,1]))
    assert c.count == 1
    assert c.base_count == 3
    assert c.mutation_count == 0
    assert c.qual == [1,1,1]
    assert c.seq == "ACA"

    c1 = ConsensusQSequence(("AAA", [1,2,3]))
    c1.add(ConsensusQSequence(("ATG", [1,1,3])))
    c = ConsensusQSequence(c1)
    assert c.count == 2
    assert c.base_count == 6
    assert c.mutation_count == 2
    assert c.qual == [1,2,0]
    assert c.seq == "AAN"

    class Empty:
        pass

    myobj = Empty()
    myobj.seq = "ANNAG"
    myobj.qual = [3,0,0,2,9]
    myobj.count = 22
    c = ConsensusQSequence(myobj)
    assert c.count == 22
    assert c.base_count == 22 * 5 
    assert c.mutation_count == 0
    assert c.qual == [3,0,0,2,9]
    assert c.seq == "ANNAG"

    c = ConsensusQSequence(ConsensusQSequence(("TGCG", [1,2,3,4], 39)))
    assert c.count == 39
    assert c.base_count == 39 * 4
    assert c.mutation_count == 0
    assert c.qual == [1,2,3,4]
    assert c.seq == "TGCG"

    c = ConsensusQSequence(None)
    assert c.count == 0
    assert c.base_count == 0
    assert c.mutation_count == 0
    assert c.qual == []
    assert c.seq == ""

    # Test initialization that unpickling protocol will use.
    # This particular initialization is not intended for manual use and its
    # input values are not (yet) extensively checked.
    c = ConsensusQSequence(((
    #   A  C  G  T  N
        0, 4, 3, 0, 0,
        5, 0, 5, 5, 0,
        0, 0, 0, 0, 4,
        0, 1, 0, 0, 4),( # quality scores
        0,40,41, 0, 0,
        1, 0, 3, 2, 0,
        0, 0, 0, 0, 0,
        0,40, 0, 0, 0), 5))
    assert c.count == 5
    assert c.base_count == 31
    assert c.mutation_count == 17
    assert c.qual == [40, 3, 0, 40]
    assert c.seq == "CGNC"

    # Test whether sequence count is correctly updated when adding two
    # ConsensusQSequence instances.
    c = ConsensusQSequence((  "ACTG", [40]*4, 5))
    c.add(ConsensusQSequence(("ACT" , [30]*3, 3)))
    c.add(ConsensusQSequence((   "G", [40])), offset = 3)
    assert c.count == 9
    assert c.base_count == 30
    assert c.mutation_count == 0
    assert c.qual == [40]*4
    assert c.seq == "ACTG"

    # Test exceptions upon initialisation with bad input values
    with pytest.raises(Exception):
        c = ConsensusQSequence(("AGCGT")) # missing qual

    with pytest.raises(Exception):
        c = ConsensusQSequence(("AGCGT", [1,3,4])) # wrong qual length

    with pytest.raises(Exception):
        c = ConsensusQSequence(("GNG", [1,1,1])) # N cannot have non-zero qual

    with pytest.raises(Exception):
        c = ConsensusQSequence(("GCC", [30,30,42])) # Q score > 41

    with pytest.raises(Exception):
        c = ConsensusQSequence(("GCC", [30,30,-1])) # Q score < 0

    with pytest.raises(Exception):
        c = ConsensusQSequence(("GQC", [30,30,30])) # non ACTGN character

    # Add a sequence
    c = ConsensusQSequence()
    c.add(ConsensusQSequence(("AAAA", [30, 30, 30, 30])))
    assert c.count == 1
    assert c.base_count == 4
    assert c.mutation_count == 0
    assert c.qual == [30]*4
    assert c.seq == "A"*4

    # Test for ambiguous (N) character and updated quality values
    c.add(ConsensusQSequence(("TAAA", [30, 29, 31, 30])))
    assert c.count == 2
    assert c.base_count == 8 
    assert c.mutation_count == 1
    assert c.qual == [0, 30, 31, 30]
    assert c.seq == "NAAA"

    # Test if most abundant character is selected and quality values are
    # updated appropriately.
    c.add(ConsensusQSequence(("TAAG", [35, 1, 1, 40])))
    assert c.count == 3
    assert c.base_count == 12 
    assert c.mutation_count == 2 
    assert c.qual == [35, 30, 31, 30]
    assert c.seq == "TAAA"

    # Test if ties on abundance are broken by quality value
    c.add(ConsensusQSequence(("TAAGG", [40]*5)))
    assert c.count == 4
    assert c.base_count == 17 
    assert c.mutation_count == 3 
    assert c.qual == [40, 40, 40, 40, 40]
    assert c.seq == "TAAGG"

    # Test adding with different offset
    c.add(ConsensusQSequence(("AAGT", [40]*4)), offset = 1)
    assert c.count == 5
    assert c.base_count == 21 
    assert c.mutation_count == 4 
    assert c.qual == [40, 40, 40, 40, 0]
    assert c.seq == "TAAGN"

    c.add(ConsensusQSequence(("CCCTAAGGA", [30]*3 + [40]*6)), offset = -3)
    assert c.count == 6
    assert c.base_count == 30 
    assert c.mutation_count == 4
    assert c.qual == [30]*3 + [40]*6
    assert c.seq == "CCCTAAGGA"

    c.add(ConsensusQSequence(("TTA", [30, 41, 1])), offset = 2)
    assert c.count == 7
    assert c.base_count == 33
    assert c.mutation_count == 5
    assert c.qual == [30, 30, 0, 41] + [40]*5
    assert c.seq == "CCNTAAGGA"

    # Test appending
    c.add(ConsensusQSequence(("GGAC", [40]*4)), offset = len(c.seq))
    assert c.count == 8
    assert c.base_count == 37
    assert c.mutation_count == 5
    assert c.qual == [30, 30, 0, 41] + [40]*9
    assert c.seq == "CCNTAAGGAGGAC"

    # Test prepending
    c.add(ConsensusQSequence(("TTT", [10]*3, 2)), offset = -3)
    assert c.count == 10
    assert c.base_count == 43
    assert c.mutation_count == 5
    assert c.qual == [10]*3 + [30, 30, 0, 41] + [40]*9
    assert c.seq == "TTTCCNTAAGGAGGAC"

    # Test adding with too large offset
    with pytest.raises(ValueError):
        c.add(ConsensusQSequence(("ACC", [40]*3)), offset = len(c.seq)+1)

    # Test adding with too small offset
    with pytest.raises(ValueError):
        c.add(ConsensusQSequence(("GGAAC", [40]*5)), offset = -6)

    # Test that failing to add did not modify anything
    assert c.count == 10
    assert c.base_count == 43
    assert c.mutation_count == 5
    assert c.qual == [10]*3 + [30, 30, 0, 41] + [40]*9
    assert c.seq == "TTTCCNTAAGGAGGAC"
    
    # Test slightly more involved combination: 
    # s2 will have offset 10, include 4 mismatches and win with 1 of them on
    # quality, and lose with another. So 2 N should be in the consensus and one
    # char changed from before. Also testing if order of adding matters.
    s1 = "ACCATCTAGAGCCGACTATCAGCAGCATCATATTACTAAGCAGCATACGACAGA"
    q1 = [10]*len(s1)
    #                      *   *                  *         *
    s2 =           "GCCGACTGTCATCAGCATCATATTACTAAGGAGCATACGATAGA"
    q2 = [10]*7 + [9] + [10]*3 + [11] + [10]*32
    s3 = "ACCATCTAGAGCCGACTATCATCAGCATCATATTACTAAGNAGCATACGANAGA"
    for t1, t2 in permutations([(s1, q1, 5), (s2, q2, 5)]):
        if t1[0] == s2:
            offset = -10
        else:
            offset = 10
        c1 = ConsensusQSequence(t1)
        c2 = ConsensusQSequence(t2)
        c1.add(c2, offset)
        assert c1.count == 10
        assert c1.base_count == 5*len(s1) + 5*len(s2)
        assert c1.mutation_count == 20
        assert c1.qual == [10]*21 + [11] + [10]*18 + [0] + [10]*9 + [0] + \
                [10]*3
        assert c1.seq == s3

    # Test pickling and unpickling
    for obj in [c, c1, c2]:
        s = pickle.dumps(obj)
        u = pickle.loads(s)
        assert not u is obj
        assert u.count == obj.count
        assert u.base_count == obj.base_count
        assert u.mutation_count == obj.mutation_count
        assert u.qual == obj.qual
        assert u.seq == obj.seq

def test_ConsensusQSequence_python():
    helperfunc_test_ConsensusQSequence(load_c = False)

def test_ConsensusQSequence_c():
    helperfunc_test_ConsensusQSequence(load_c = True)

def test_Clone():
    import pickle
    from itertools import permutations
    from rtcr.clone import Clone

    clone1 = Clone(("ACTG", [30]*4, 5))
    assert clone1.count == 5
    assert clone1.base_count == 20
    assert clone1.mutation_count == 0
    assert clone1.orig_count == 5
    assert clone1.qual == [30]*4
    assert len(clone1) == 4
    assert clone1.seq == "ACTG"
    assert str(clone1) == clone1.seq
    assert clone1 == clone1

    clone_dict = {clone1 : clone1}
    
    # create two more clones that differ only by quality or abundance from
    # clone1
    clone1a = Clone(("ACTG", [30]*4, 10))
    clone1b = Clone(("ACTG", [31]*4, 5))
    assert clone1a == clone1 and not clone1a != clone1
    assert clone1b == clone1 and not clone1b != clone1
    assert clone1 in clone_dict
    assert clone1a in clone_dict
    assert clone1b in clone_dict
    
    #sanity check the test itself
    assert clone1a.seq == clone1.seq
    assert clone1a.qual == clone1.qual
    assert clone1a.count != clone1.count
    assert clone1b.seq == clone1.seq
    assert clone1b.qual != clone1.qual
    assert clone1b.count == clone1.count

    clone1c = Clone(("ACTG", [30]*4, 5))
    assert clone1c == clone1

    clone2 = Clone(("CCTG", [31]*4, 3))
    # Test order of adding clones does not matter
    for a, b in permutations((clone1, clone2)):
        cc1 = a.add(b)
        # Check that clone1 did not change
        assert clone1.count == 5
        assert clone1.base_count == 20
        assert clone1.mutation_count == 0
        assert clone1.orig_count == 5
        assert clone1.qual == [30]*4
        assert len(clone1) == 4
        assert clone1.seq == "ACTG"
        # Check that cc1 is the combined clone
        assert cc1.count == 8
        assert cc1.base_count == 32
        assert cc1.mutation_count == 3
        assert cc1.orig_count == 5
        assert cc1.qual == [30, 31, 31, 31]
        assert cc1.seq == "ACTG"

    clone3 = Clone(("CCGG", [10, 40, 31, 40], 3))
    for a, b in permutations((cc1, clone3)):
        cc2 = a.add(b)
        assert cc2.count == 11
        assert cc2.base_count == 44
        assert cc2.mutation_count == 8 
        assert cc2.orig_count == 5
        assert cc2.qual == [31, 40, 31, 40]
        assert cc2.seq == "CCTG"

    cc3 = clone2.add(clone3)
    assert cc3.count ==  6
    assert cc3.base_count == 24
    assert cc3.mutation_count == 3
    assert cc3.orig_count == 3
    assert cc3.qual == [31, 40, 0, 40]
    assert cc3.seq == "CCNG"

    # Test adding clones with different reference positions
    clone4 = Clone(("ACCT", [30]*4), refpos = 2)
    clone5 = Clone((  "CGAAG", [30]*5), refpos = 0)
    cc4 = clone4.add(clone5)
    assert cc4.count == 2
    assert cc4.base_count == 9
    assert cc4.mutation_count == 1
    assert cc4.orig_count == 1
    assert cc4.qual == [30, 30, 30, 0] + [30]*3
    assert cc4.seq == "ACCNAAG"

    # Test instantiating a clone from another clone
    cc5 = Clone(cc3)
    assert cc5.count == cc3.count
    assert cc5.base_count == cc3.base_count
    assert cc5.mutation_count == cc3.mutation_count
    assert cc5.orig_count == cc3.orig_count
    assert cc5.qual == cc3.qual
    assert cc5.seq == cc3.seq

    # Test pickling
    for clone in [clone1, clone2, clone3, cc1, cc2, cc3]:
        s = pickle.dumps(clone)
        u = pickle.loads(s)
        assert not u is clone
        assert u.count == clone.count
        assert u.base_count == clone.base_count
        assert u.mutation_count == clone.mutation_count
        assert u.orig_count == clone.orig_count
        assert u.qual == clone.qual
        assert u.seq == clone.seq

def test_CloneSet():
    import pickle
    from itertools import izip
    from rtcr.clone import Clone, CloneSet

    clone1 = Clone(("ACTG", [30]*4, 5))
    cloneset = CloneSet()
    assert len(cloneset) == 0
    assert cloneset.base_count == 0
    assert cloneset.mutation_count == 0
    assert cloneset.sequence_count == 0
    
    cloneset.add(clone1)
    assert len(cloneset) == 1
    assert clone1 in cloneset
    assert cloneset.base_count == 20
    assert cloneset.mutation_count == 0
    assert cloneset.sequence_count == 5 

    with pytest.raises(Exception) as excinfo:
        cloneset.add("hello world")

    assert cloneset.count == 1

    # Test the adding of other clone objects that should not end up in the set
    # because they have the same sequence as an existing clone.
    clone1a = Clone(("ACTG", [30]*4, 10))
    clone1b = Clone(("ACTG", [31]*4, 5))
    clone1c = Clone(("ACTG", [30]*4, 5))
    assert clone1a.seq == clone1.seq
    assert clone1a.qual == clone1.qual
    assert clone1a.count != clone1.count

    assert clone1b.seq == clone1.seq
    assert clone1b.qual != clone1.qual
    assert clone1b.count == clone1.count
    
    assert clone1c.seq == clone1.seq
    assert clone1c.qual == clone1.qual
    assert clone1c.count == clone1.count

    with pytest.raises(KeyError):
        cloneset.add(clone1a)
        
    # Although technically the clone1a/b/c clones are not identical to clone1,
    # a clone's identity is considered fully determined by its sequence. Also,
    # a clone is immutable.
    assert clone1a in cloneset
    assert clone1b in cloneset
    assert clone1c in cloneset

    # Test merging of clones
    cloneset = CloneSet()
    clone1 = Clone(("AAAA", [30]*4, 100))
    clone2 = Clone(("ATTA", [30]*4, 10))
    clone3 = Clone(("ACAA", [30]*4, 10))
    clone4 = Clone(("CGCG", [25, 25, 25, 35], 1))
    clone5 = Clone(("CCGG", [30]*4, 100))

    cloneset.add(clone1)
    cloneset.add(clone2)
    cloneset.add(clone3)
    cloneset.add(clone4)
    cloneset.add(clone5)

    assert len(cloneset) == 5
    diff, clone = cloneset.merge(clone2, clone3)
    #print "clone2 = %s ; %s"%(clone2, clone2.mutation_count)
    #print "clone3 = %s ; %s"%(clone3, clone3.mutation_count)
    #print "merged = %s ; %s"%(clone, clone.mutation_count)
    assert diff == 20
    assert clone.count == 20
    assert clone.orig_count == 10
    assert clone.mutation_count == 20
    assert clone.seq == "ANNA"
    assert clone.qual == [30, 0, 0, 30]
    assert len(cloneset) == 4
    assert cloneset.base_count == 884 
    assert cloneset.mutation_count == 20
    assert cloneset.sequence_count == 221 

    diff, clone = cloneset.merge(clone1, clone)
    assert diff == 10
    assert clone.count == 120
    assert clone.orig_count == 100
    assert clone.mutation_count == 30
    assert clone.seq == "AAAA"
    assert clone.qual == [30]*4
    assert len(cloneset) == 3
    assert cloneset.base_count == 884 
    assert cloneset.mutation_count == 30
    assert cloneset.sequence_count == 221

    diff, clone = cloneset.merge(clone5, clone4)
    assert diff == 2 
    assert clone.count == 101
    assert clone.orig_count == 100
    assert clone.mutation_count == 2
    assert clone.seq == "CCGG"
    assert clone.qual == [30, 30, 30, 35]
    assert len(cloneset) == 2
    assert cloneset.base_count == 884 
    assert cloneset.mutation_count == 32
    assert cloneset.sequence_count == 221

    # Test merging clones, updating their orig_count fields
    cs = CloneSet()
    c1 = Clone(("ACT", [40]*3, 5), orig_count = 2)
    assert c1.orig_count == 2
    c2 = Clone(("ACT", [30]*3, 7), orig_count = 4)
    assert c2.orig_count == 4
    cs.add(c1)
    assert cs.sequence_count == 5
    with pytest.raises(KeyError):
        cs.add(c2)
    assert cs.sequence_count == 5
    cs.add(c2, merge = True)
    c = next(iter(cs))
    assert c.count == 12
    assert c.base_count == 36
    assert c.mutation_count == 0
    assert c.orig_count == 6

    # TODO: test adding clonesets together

    # Test pickling
    s = pickle.dumps(cloneset)
    u = pickle.loads(s)
    assert not u is cloneset
    assert len(u) == len(cloneset)
    assert u.base_count == cloneset.base_count
    assert u.mutation_count == cloneset.mutation_count
    assert u.sequence_count == cloneset.sequence_count
    for uc, c in izip(u, cloneset):
        assert not uc is c
        assert uc.count == c.count
        assert uc.orig_count == c.orig_count
        assert uc.mutation_count == c.mutation_count
        assert uc.seq == c.seq
        assert uc.qual == c.qual

def helperfunc_test_simple_clone_merging(gen_clone):
    from rtcr.clone import CloneSet
    c0 = gen_clone(("A"*45, [40]*45, 63))
    c1 = gen_clone(("A"*44 + "T", [40]*45, 29))
    c2 = gen_clone(("A"*43 + "TT", [40]*45, 6))
    c3 = gen_clone(("A"*42 + "TTT", [40]*45, 2))

    # Test if order of adding clones matters (it should not!), first on the
    # Clone object level.
    # First ordering: 2->1, 1->0, 3->0
    c_2_1 = c1.add(c2)
    assert c_2_1.count == 35
    assert c_2_1.base_count == 35*45
    assert c_2_1.mutation_count == 6
    assert c_2_1.seq == c1.seq
    c_1_0 = c0.add(c_2_1)
    assert c_1_0.count == 98
    assert c_1_0.base_count == 98*45
    assert c_1_0.mutation_count == 41
    assert c_1_0.seq == c0.seq
    c_3_0 = c_1_0.add(c3)
    assert c_3_0.count == 100
    assert c_3_0.base_count == 4500
    assert c_3_0.mutation_count == 47 
    assert c_3_0.seq == c0.seq
    # Second ordering: 1->0, 3->2, 2->0
    c_1_0 = c0.add(c1)
    assert c_1_0.count == 92
    assert c_1_0.base_count == 92*45
    assert c_1_0.mutation_count == 29
    assert c_1_0.seq == c0.seq
    c_3_2 = c2.add(c3)
    assert c_3_2.count == 8
    assert c_3_2.base_count == 8*45
    assert c_3_2.mutation_count == 2
    assert c_3_2.seq == c2.seq
    c_2_0 = c_1_0.add(c_3_2)
    assert c_2_0.count == 100
    assert c_2_0.base_count == 4500
    assert c_2_0.mutation_count == 47
    assert c_2_0.seq == c0.seq

    #Test the above merge processes on CloneSet
    cs = CloneSet([c0, c1, c2, c3])
    assert cs.count == 4
    assert cs.sequence_count == 100
    assert cs.base_count == 4500
    assert cs.mutation_count == 0
    # First ordering: 2->1, 1->0, 3->0
    mc_diff, c_2_1 = cs.premerge(c1, c2)
    assert mc_diff == 6
    assert c_2_1.count == 35
    assert c_2_1.base_count == 35*45
    assert c_2_1.mutation_count == 6
    assert c_2_1.seq == c1.seq
    # check this premerge did not affect the cloneset yet
    assert cs.count == 4
    assert cs.sequence_count == 100
    assert cs.base_count == 4500
    assert cs.mutation_count == 0
    # Perform the merge
    cs.enact_merge(c1, c2, c_2_1)
    assert cs.count == 3
    assert cs.sequence_count == 100
    assert cs.base_count == 4500
    assert cs.mutation_count == 6
    # NOTE: using c1 and not c_2_1 to check if the cloneset correctly uses its
    # own clone rather than the one provided!
    mc_diff, c_1_0 = cs.premerge(c1, c0)
    assert mc_diff == 35
    assert c_1_0.count == 98
    assert c_1_0.base_count == 98*45
    assert c_1_0.mutation_count == 41
    assert c_1_0.seq == c0.seq
    cs.enact_merge(c1, c0, c_1_0)
    assert cs.count == 2
    assert cs.sequence_count == 100
    assert cs.base_count == 4500
    assert cs.mutation_count == 41
    mc_diff, c_3_0 = cs.premerge(c3, c0)
    assert mc_diff == 6
    assert c_3_0.count == 100
    assert c_3_0.base_count == 4500
    assert c_3_0.mutation_count == 47 
    assert c_3_0.seq == c0.seq
    cs.enact_merge(c3, c0, c_3_0)
    assert cs.count == 1 
    assert cs.sequence_count == 100
    assert cs.base_count == 4500
    assert cs.mutation_count == 47
    clone = next(iter(cs))
    assert clone.count == 100
    assert clone.base_count == 4500
    assert clone.mutation_count == 47
    assert clone.seq == c0.seq

    # Second ordering: 1->0, 3->2, 2->0
    cs = CloneSet([c0, c1, c2, c3])
    assert cs.count == 4
    assert cs.sequence_count == 100
    assert cs.base_count == 4500
    assert cs.mutation_count == 0
    mc_diff, c_1_0 = cs.premerge(c1, c0)
    assert mc_diff == 29
    assert c_1_0.count == 92
    assert c_1_0.base_count == 92*45
    assert c_1_0.mutation_count == 29
    assert c_1_0.seq == c0.seq
    cs.enact_merge(c1, c0, c_1_0)
    assert cs.count == 3
    assert cs.sequence_count == 100
    assert cs.base_count == 4500
    assert cs.mutation_count == 29
    mc_diff, c_3_2 = cs.premerge(c3, c2)
    assert mc_diff == 2
    assert c_3_2.count == 8
    assert c_3_2.base_count == 8*45
    assert c_3_2.mutation_count == 2
    assert c_3_2.seq == c2.seq
    cs.enact_merge(c3, c2, c_3_2)
    assert cs.count == 2 
    assert cs.sequence_count == 100
    assert cs.base_count == 4500
    assert cs.mutation_count == 31
    mc_diff, c_2_0 = cs.premerge(c2, c0)
    assert mc_diff == 16
    assert c_2_0.count == 100
    assert c_2_0.base_count == 4500
    assert c_2_0.mutation_count == 47
    assert c_2_0.seq == c0.seq
    cs.enact_merge(c2, c0, c_2_0)
    assert cs.count == 1 
    assert cs.sequence_count == 100
    assert cs.base_count == 4500
    assert cs.mutation_count == 47

def test_simple_Clone_merging():
    from rtcr.clone import Clone
    def gen_clone(*args, **kwargs):
        return Clone(*args, **kwargs)
    helperfunc_test_simple_clone_merging(gen_clone)

def test_simple_AnnotatedClone_merging(vjref):
    from rtcr.seq import IntervalAnnotation
    from rtcr.clone import AnnotatedClone
    vid = "TRAV1*01"
    jid = "TRAJ1*01"
    v_annot = IntervalAnnotation(
        allele = vjref[vid],
        start = 0,
        end = 10,
        score = 20)

    j_annot = IntervalAnnotation(
        allele = vjref[jid],
        start = 15,
        end = 25,
        score = 20)
    def gen_clone(*args, **kwargs):
        return AnnotatedClone(v_annot, None, j_annot, None, *args, **kwargs)
    helperfunc_test_simple_clone_merging(gen_clone)

def test_Trie_python_specific():
    from rtcr.trie.trie import Trie

    t = Trie()
    for s in ["{:08b}".format(i) for i in xrange(256)]:
        t[s] = 0 # "00000000" to "11111111"
 
    # Test selecting certain nodes
    for s in ["{:04b}".format(i) for i in xrange(16)]:
        t[s] = 1 # "0000" to "1111"
    assert len(list(t.pairs_ext(1, lambda node:node.value==1))) == (16*4)/2

    ######################################
    # Test get_nearest_variants function #
    ######################################
    t = Trie()
    for s in ["AAAAAAAAAA",
            "AAAAAAAATA",
            "AAAAATAATA",
            "AAAAATAAAA",
            "TAAAAAAAAA",
            "ATAAAAAAAA"]:
        t[s] = 0
    assert set([(hd, s, value) \
            for hd, s, value in t.get_nearest_variants("A"*10)]) == set([
        (1, "AAAAAAAATA", 0),
        (1, "AAAAATAAAA", 0),
        (1, "TAAAAAAAAA", 0),
        (1, "ATAAAAAAAA", 0)])

    assert set([(hd, s, value) \
            for hd, s, value in t.get_nearest_variants("AAAAATAATA")]) == \
            set([(1, "AAAAAAAATA", 0),
                (1, "AAAAATAAAA", 0)])

    assert set([(hd, s, value) \
            for hd, s, value in t.get_nearest_variants("A"*11)]) == set()

    assert set([(hd, s, value) \
            for hd, s, value in t.get_nearest_variants("AAAAATACTC")]) == \
            set([(2, "AAAAATAATA", 0)])
    
    # Test maxhd
    assert set([(hd, s, value) \
            for hd, s, value in t.get_nearest_variants("AAAAATACTC",
                maxhd = 1)]) == set()

def eqp(pairs, x):
    """test if pairs and x are equal."""
    if len(pairs) != len(x):
        return False
    for pair in pairs:
        if not (pair in x or pair[::-1] in x):
            return False
    return True

def test_SearchableCloneSet():
    from itertools import izip
    from collections import deque
    from rtcr.clone import Clone, CloneSet, SearchableCloneSet

    clone = Clone(("ACTG", [30]*4, 5))
    cloneset = CloneSet()
    assert len(cloneset) == 0

    scls = SearchableCloneSet(cloneset)
    assert len(scls) == 0

    scls.add(clone)
    assert scls._trie[clone.seq] == [clone] 
    assert len(scls) == 1
    assert len(cloneset) == 0

    scls.remove(clone)
    with pytest.raises(KeyError):
        scls._trie[clone.seq]
    assert len(scls) == 0

    gen_clone = lambda seq, count: Clone((seq, [40]*len(seq), count))

    c0 = gen_clone("AAAAAAAAAAAAAAAAAAAA", 67)
    c1 = gen_clone("AAAAAATAAAAAAAAAAAAA", 5)
    c2 = gen_clone("AAATAAAAAAAAAATAAAAA", 2)
    c3 = gen_clone("AATAATAAAATAATAATAAA", 1)
    cloneset = CloneSet((c0, c1, c2, c3))
    # Test searching for clones
    scls = SearchableCloneSet(cloneset)

    assert set([(hd, clone.seq, clone.count) for hd, clone in \
            scls.neighbors("AAAAAAAAAAAAAAAAAAAA", maxhd = 2)]) == \
            set([(1, "AAAAAATAAAAAAAAAAAAA", 5),
            (2, "AAATAAAAAAAAAATAAAAA", 2)])

    assert set([(hd, clone.seq, clone.count) for hd, clone in \
            scls.neighbors(c0, maxhd = 2)]) == \
            set([(1, "AAAAAATAAAAAAAAAAAAA", 5),
            (2, "AAATAAAAAAAAAATAAAAA", 2)])

    # Test nearest_neighbors
    #assert set([(hd, clone.seq, clone.count) for hd, clone in \
    #        scls.nearest_neighbors(c0)]) == \
    #        set([(1, "AAAAAA*AAAAAAAAAAAAA", 5)])

    #assert set([(hd, clone.seq, clone.count) for hd, clone in \
    #        scls.nearest_neighbors(c3, maxhd = 5)]) == \
    #        set([(5, c0.seq, 67)])

    #assert set([(hd, clone.seq, clone.count) for hd, clone in \
    #        scls.nearest_neighbors(c3, maxhd = 4)]) == set()

    # Test combination of merging and searching for clones
    assert len(scls) == 4
    assert scls.sequence_count == 75
    assert set([(hd, clone.seq, clone.count) for hd, clone in \
            scls.neighbors(c0.seq, maxhd = 5)]) == \
            set([(1, "AAAAAATAAAAAAAAAAAAA", 5),
            (2, "AAATAAAAAAAAAATAAAAA", 2),
            (5, "AATAATAAAATAATAATAAA", 1)])

    diff, clone = scls.merge(c0, c3)
    assert diff == 5
    assert clone.count == 68
    assert clone.orig_count == 67
    assert clone.mutation_count == 5
    assert clone.seq == c0.seq
    assert clone.qual == [40]*len(c0.seq)
    assert len(scls) == 3
    assert scls.base_count == 20*75 
    assert scls.mutation_count == 5 
    assert scls.sequence_count == 75
    assert set([(hd, clone.seq, clone.count) for hd, clone in \
            scls.neighbors(c0.seq, maxhd = 5)]) == \
            set([(1, "AAAAAATAAAAAAAAAAAAA", 5),
            (2, "AAATAAAAAAAAAATAAAAA", 2)])
    diff, clone = scls.merge(c0, c1)
    assert set([(hd, clone.seq, clone.count) for hd, clone in \
            scls.neighbors(c0.seq, maxhd = 5)]) == \
            set([(2, "AAATAAAAAAAAAATAAAAA", 2)])
    diff, clone = scls.merge(c0, c2)
    assert list(scls.neighbors(c0.seq, maxhd = 5)) == []

    seq2clone = lambda seq: Clone((seq, [40]*len(seq), 1))
    # Test searching for clone pairs
    scls2 = SearchableCloneSet(CloneSet())
    scls2.add(seq2clone("AAAA"))
    scls2.add(seq2clone("AAAT"))
    scls2.add(seq2clone("TAAT"))
    scls2.add(seq2clone("TATA"))

    get_pairs = lambda scls, seqlen, maxhd: [(clone1.seq, clone2.seq) \
            for (hd, clone1, clone2) in scls.pairs(seqlen = seqlen,
                maxhd = maxhd)]
    assert eqp(get_pairs(scls2, 4, 1), [("AAAA", "AAAT"),
            ("AAAT", "TAAT")])
    for hd, clone1, clone2 in scls2.pairs(4, 2):
        assert hd == sum([ch1 != ch2 \
                for ch1, ch2 in izip(clone1.seq, clone2.seq)])

def helperfunc_test_get_offset(load_c):
    if load_c:
        from cseq import ConsensusQSequence, get_offset
    else:
        from rtcr.seq.seq import ConsensusQSequence, get_offset

    s1 =   "ACTG"
    s2 = "AAACTG"
    offset, score = get_offset(s1, s2, k = 4, moff = 4)
    assert offset == 2
    assert score == 4

    # Test scoring function
    # (mis)match with "N" is 0, match is 1, mismatch is -1
    s1 =   "ANTG"
    s2 = "AAACTG"
    offset, score = get_offset(s1, s2, k = 4, moff = 4)
    assert offset == 2
    assert score == 3
    offset, score = get_offset("AGTG", s2, k = 4, moff = 4)
    assert offset == 2
    assert score == 2
    offset, score = get_offset("ANTG", "AAANTG", k = 4, moff = 4)
    assert offset == 2
    assert score == 3

    # get_offset method chooses a k-mer in the middle of the longer sequence,
    # (if tied, chooses s1).
    #                       __________ (location of 10-mer in s1)
    s1 =          "TGTGCGTGACAGATAAAGCGCTGAGCTG"
    s2 = "ACCAAGGCTTGTGCGTGACAGATAAAGC"
    offset, score = get_offset(s1, s2, k = 10, moff = 12)
    assert offset == 9
    assert score == 10
    s1 = "ACTGATCATCCGAAGCCGGCTGTGGCGCATCGACGGC"
    s2 =       "CATCCGAAGCCGGCTGTGGCGCATCGACGGCGCGTG"
    offset, score = get_offset(s1, s2, k = 10, moff = 12)
    assert offset == -6
    assert score == 10

    # This test revealed a mistake in one of the get_offset C implementations
    # where s1 and s2 were swapped, but not their lengths.
    s1 = "AGGCTTGTGCGTGACAGATAAAGCGCTGAGCTGAGCTCTG"
    s2 =      "TGTGCGTGACAGATAAAGCGCTGAGCTG"
    assert get_offset(s1, s2, k = 10, moff = 5) == (-5, 10)

    #TODO: test edge cases
    #TODO: test when None value should be returned

def test_get_offset_c():
    helperfunc_test_get_offset(load_c = True)

def test_get_offset_python():
    helperfunc_test_get_offset(load_c = False)

def test_umi_group_ec():
    from itertools import permutations
    from rtcr.seq import QSequence
    from rtcr.umi_group_ec import run_umi_group_ec
    
    r1 = QSequence("UMI:R1:A:I", "AGGCTTGTGCGTGACAGATAAAGCGCTGAGCTGAGCTCTG",
            [40]*40)
    r2 = QSequence("UMI:R1:A:I",      "TGTGCGTGACAGATAAAGCGCTGAGCTG", [40]*28)

    class DummyWriter(object):
        def __init__(self):
            self.records = []
        def __add__(self, record):
            self.records += [record]
            return self
        def __len__(self):
            return len(self.records)

    for records in [(r1, r2), (r2, r1)]:
        discarded = DummyWriter()
        res = run_umi_group_ec(records, k = 10, moff = 5,
                min_score4offset = 0, min_score4merge = 0,
                discarded = discarded)
        assert len(discarded) == 0
        assert len(res) == 1
        assert len(res['R1:A']) == 1
        rec = res['R1:A'][0]
        assert rec.seq == r1.seq
        assert rec.qual == [40]*40
        assert rec.count == 2
        assert rec.base_count == 68
        assert rec.mutation_count == 0

    # r3 has 4 extra characters prepended to start of r1
    r3 = QSequence('UMI:R1:A:I', "ACCAAGGCTTGTGCGTGACAGATAAAGCGCTGAGCT",
            [40]*36)
    for records in permutations((r1, r2, r3)):
        discarded = DummyWriter()
        res = run_umi_group_ec(records, k = 15, moff = 10,
                min_score4offset = 0, min_score4merge = 0,
                discarded = discarded)
        assert len(discarded) == 0
        assert len(res) == 1
        assert len(res['R1:A']) == 1
        rec = res['R1:A'][0]
        assert rec.seq == "ACCA" + r1.seq
        assert rec.qual == [40]*44
        assert rec.count == 3 
        assert rec.base_count == 104 
        assert rec.mutation_count == 0

    discarded = DummyWriter()
    res = run_umi_group_ec((r1, r2, r3), k = 50, moff = 5,
            min_score4offset = 0, min_score4merge = 0,
            discarded = discarded)
    assert len(discarded) == 3 
    assert len(res) == 0

    r4 = QSequence('UMI:R1:T:I',  'ACCCCAAAAAGGTGGTG', [40]*17)
    r5 = QSequence('UMI:R1:T:I', 'AACCCCAAAAAGGGGG', [40]*16)
    for records in permutations((r1, r2, r3, r4, r5)):
        discarded = DummyWriter()
        res = run_umi_group_ec(records, k = 15, moff = 10,
                min_score4offset = 0, min_score4merge = 0,
                discarded = discarded)
        assert len(discarded) == 0
        assert len(res) == 2
        assert len(res['R1:A']) == 1
        assert len(res['R1:T']) == 1
        rec_a = res['R1:A'][0]
        rec_b = res['R1:T'][0]
        assert rec_a.seq == "ACCA" + r1.seq
        assert rec_a.qual == [40]*44
        assert rec_a.count == 3
        assert rec_a.base_count == 104 
        assert rec_a.mutation_count == 0
        assert rec_b.seq == 'AACCCCAAAAAGGNGGTG'
        assert rec_b.qual == [40]*13 + [0] + [40]*4
        assert rec_b.count == 2
        assert rec_b.base_count == 33 
        assert rec_b.mutation_count == 1

    # TODO: more extensive tests with larger datasets
    # TODO: test min_score4offset parameter
    # TODO: test min_score4merge parameter
    # TODO: test UMI qualitifer (e.g. 'R1')

def test_rtcr_heap():
    from itertools import permutations
    from rtcr.heap import Heap
    h = Heap()

    push = lambda iterable:[h.push(it) for it in iterable]
    pop = lambda:[h.pop() for i in xrange(len(h))]
    popjoin = lambda:"".join(pop())

    # Test if numbers are popped in the correct order 
    assert len(h) == 0
    push(range(10))
    assert len(h) == 10
    assert pop() == range(10)
    push(range(10)[::-1])
    assert len(h) == 10
    assert pop() == range(10)

    # Test deletion of specific items
    assert len(h) == 0
    push("ACGDFEB")
    assert popjoin() == "ABCDEFG"
    assert len(h) == 0
    push("ACGDFEB")
    h.remove_item("D")
    assert popjoin() == "ABCEFG"
    push("ACGDFEB")
    h.remove_item("E")
    assert popjoin() == "ABCDFG"
    push("ACGDFEB")
    h.remove_item("B")
    assert popjoin() == "ACDEFG"
    push("ACGDFEB")
    h.remove_item("G")
    assert popjoin() == "ABCDEF"

    # Test has_item method
    assert len(h) == 0
    push("ACGDFEB")
    for ch in "ABCDEFG":
        assert h.has_item(ch)

    # Test peeking ahead
    assert len(h) == 7
    assert h.peek() == "A"
    while len(h) > 0:
        assert h.peek() == h.pop()

    # Test adding the same item twice fails
    with pytest.raises(Exception) as excinfo:
        assert len(h) == 0
        h.push("A")
        h.push("A")
    assert excinfo.value.message == "Item already on the heap."

    assert len(h) == 1
    assert h.pop() == "A"

    # Test popping empty heap fails
    with pytest.raises(Exception) as excinfo:
        assert len(h) == 0 
        h.pop()
    assert excinfo.value.message == "Cannot pop from empty heap."

    # Test peeking empty heap fails
    with pytest.raises(Exception) as excinfo:
        assert len(h) == 0
        h.peek()
    assert excinfo.value.message == "Cannot peek empty heap."

def helperfunc_test_cigar_intervals(load_c):
    if load_c:
        from cseq import cigar_intervals
    else:
        from rtcr.seq.seq import cigar_intervals

    # Interval tests without indels and clipping
    assert cigar_intervals("10M", 0) == (0, 10, 0, 10)
    assert cigar_intervals("10M", 5) == (0, 10, 5, 15)
    assert cigar_intervals("5M5=", 5) == (0, 10, 5, 15)
    assert cigar_intervals("5M3X2=", 5) == (0, 10, 5, 15)

    # Interval test with soft clipping
    assert cigar_intervals("3S10M5S", 5) == (3, 13, 5, 15)

    # Interval test with deletion
    assert cigar_intervals("5M1D4M", 5) == (0, 9, 5, 15)

    # Interval test with insertion
    assert cigar_intervals("5M1I5M", 5) == (0, 11, 5, 15)

    # Interval test with soft clipping and indels
    assert cigar_intervals("5S10M1D5M2I10S", 4) == (5, 22, 4, 20)

def helperfunc_test_cigar_rpos2qpos(load_c):
    if load_c:
        from cseq import cigar_rpos2qpos
    else:
        from rtcr.seq.seq import cigar_rpos2qpos

    # Test rpos2qpos without indels
    for i in xrange(5):
        assert cigar_rpos2qpos("3S5M", 0, i) == 3 + i
        assert cigar_rpos2qpos("5M", 0, i) == i
        assert cigar_rpos2qpos("5M", 10, 10 + i) == i

    # Test positions outside of alignment result in None 
    assert cigar_rpos2qpos("3S5M", 0, -1) is None
    assert cigar_rpos2qpos("3S5M", 10, 9) is None
    assert cigar_rpos2qpos("3S5M", 0, 5) is None

    # Test rpos2qpos with indels
    cigar = "3S3M3D8M3I1M5S"
    for i in xrange(3): # the 3M in the cigar
        assert cigar_rpos2qpos(cigar, 0, i) == 3 + i
    for i in xrange(3): # the 3D in the cigar
        assert cigar_rpos2qpos(cigar, 0, 3 + i) is None
    for i in xrange(8): # the 8M in the cigar
        assert cigar_rpos2qpos(cigar, 0, 6 + i) == 6 + i
    assert cigar_rpos2qpos(cigar, 0, 14) == 17 # the 1M at the end
    assert cigar_rpos2qpos(cigar, 0, 15) is None

def test_cigar_intervals_c():
    helperfunc_test_cigar_intervals(load_c = True)

def test_cigar_intervals_python():
    helperfunc_test_cigar_intervals(load_c = False)

def test_cigar_rpos2qpos_c():
    helperfunc_test_cigar_rpos2qpos(load_c = True)

def test_cigar_rpos2qpos_python():
    helperfunc_test_cigar_rpos2qpos(load_c = False) 

@pytest.fixture
def vjref():
    from itertools import product
    from rtcr.allele import Allele, AlleleContainer

    gen_allele = lambda name,seq,refpos:Allele("TestSpecies", name, "F",
            refpos, str(hash(name)), seq)
    ref = AlleleContainer()

    # Generate 25 V and 25 J sequences, each 50 bases long
    for i, p in enumerate(product(*["ACTG"]*6)):
        p = "".join(p)
        if i > 50:
            break
        if i <= 25:
            ref.add(gen_allele("TRAV%s*0%s"%((i // 5) + 1, i % 5),
                "A"*40 + p + "A"*4, 40))
        if i > 25:
            ref.add(gen_allele("TRAJ%s*0%s"%(((i-25) // 5) + 1, i % 5),
                "C"*4 + p + "C"*40, 10))
    return ref

def randseq(seqlen):
    return "".join(random.choice("ACTG") for _ in xrange(seqlen))

def __test_ConsensusIntervalAnnotation(vjref):
    from rtcr.seq import ConsensusIntervalAnnotation, IntervalAnnotation
    from rtcr.allele import Allele
    from itertools import permutations

    cia = ConsensusIntervalAnnotation()
    assert cia.allele is None
    assert cia.start == 0
    assert cia.end == 0
    assert cia.score == 0
    
    # Interval annotation allows allele to be None
    IntervalAnnotation(None, 0, 0, 0)

    with pytest.raises(ValueError):
        IntervalAnnotation(None, -1, 0, 0)

    with pytest.raises(ValueError):
        IntervalAnnotation(None, 0, -1, 0)

    with pytest.raises(ValueError):
        IntervalAnnotation(None, 0, 0, -1)

    # end is not allowed to be smaller than start
    with pytest.raises(ValueError):
        IntervalAnnotation(None, 5, 3, 0)

    vid = "TRAV5*01"
    IA1 = IntervalAnnotation(vjref[vid], 5, 15, 20)
    IA2 = IntervalAnnotation(vjref[vid], 1, 2, 21)

    # Test initialization with IntervalAnnotation
    cia = ConsensusIntervalAnnotation(IA1)
    assert cia.allele == IA1.allele
    assert cia.start == IA1.start
    assert cia.end == IA1.end
    assert cia.score == IA1.score

    # Test initialization with ConsensusIntervalAnnotation
    cia1 = ConsensusIntervalAnnotation(IA1)
    cia2 = ConsensusIntervalAnnotation(cia)
    assert cia2.allele == cia1.allele
    assert cia2.start == cia1.start
    assert cia2.end == cia1.end
    assert cia2.score == cia1.score

    # Test combination of initialization and addition
    for ia1, ia2 in permutations((IA1, IA2)):
        cia1 = ConsensusIntervalAnnotation(ia1).add(ia2)
        cia2 = ConsensusIntervalAnnotation()
        cia2.add(ia1)
        cia2.add(ia2)
        assert cia2.allele == cia1.allele
        assert cia2.start == cia1.start
        assert cia2.end == cia1.end
        assert cia2.score == cia1.score

    # Test interval gets updated based on score
    for ia1, ia2 in [(IA1, IA2), (IA2, IA1)]:
        cia = ConsensusIntervalAnnotation()
        cia.add(ia1)
        assert cia.allele.name == ia1.allele.name 
        assert cia.start == ia1.start 
        assert cia.end == ia1.end
        assert cia.score == ia1.score
        cia.add(ia2)
        assert cia.allele.name == vid
        assert cia.start == IA2.start
        assert cia.end == IA2.end

    # Test if frequency trumps score for interval
    IA3 = IntervalAnnotation(vjref[vid], 3, 5, 19)
    for p in permutations(((IA1, 1), (IA2, 1), (IA3, 2))):
        cia = ConsensusIntervalAnnotation()
        for ia, count in p:
            cia.add(ia, count)
        assert cia.allele.name == vid
        assert cia.start == IA3.start
        assert cia.end == IA3.end

    vid1 = "TRAV5*01"
    vid2 = "TRAV5*02"
    vid3 = "TRAV5*03"
    # Test if allele gets updated based on score
    IA4 = IntervalAnnotation(vjref[vid1], 5, 15, 20)
    IA5 = IntervalAnnotation(vjref[vid2], 5, 15, 21)
    for ia1, ia2 in permutations((IA4, IA5)):
        cia = ConsensusIntervalAnnotation()
        cia.add(ia1)
        assert cia.allele.name == ia1.allele.name 
        assert cia.start == ia1.start 
        assert cia.end == ia1.end
        assert cia.score == ia1.score
        cia.add(ia2)
        assert cia.allele.name == IA5.allele.name 
        assert cia.start == IA5.start 
        assert cia.end == IA5.end

    # Test if frequency trumps score for allele
    IA6 = IntervalAnnotation(vjref[vid3], 5, 15, 21)
    for p in permutations(((IA4, 1), (IA5, 1), (IA6, 2))):
        cia = ConsensusIntervalAnnotation()
        for ia, count in p:
            cia.add(ia, count)
        assert cia.allele.name == IA6.allele.name 
        assert cia.start == IA6.start
        assert cia.end == IA6.end

    # Test if allele trumps interval
    IA7 = IntervalAnnotation(vjref[vid1], 5, 15, 20)
    IA8 = IntervalAnnotation(vjref[vid2], 5, 15, 20)
    IA9 = IntervalAnnotation(vjref[vid3], 0, 0, 20)
    for p in permutations(((IA7, 1), (IA8, 1), (IA9, 2))):
        cia = ConsensusIntervalAnnotation()
        for ia, count in p:
            cia.add(ia, count)
        assert cia.allele.name == IA9.allele.name 
        assert cia.start == IA9.start
        assert cia.end == IA9.end

    # Test we cannot add allele from another gene or another region
    cia = ConsensusIntervalAnnotation()
    trav2 = vjref["TRAV2*01"]
    trbv2 = Allele(species = trav2.species,
            name = "TRBV2*01", functionality = trav2.functionality,
            refpos = trav2.refpos, acc_nr = str(hash("TRBV2*01")),
            seq = trav2.seq)
    cia.add(IntervalAnnotation(trav2, 0, 0, 0))
    cia.add(IntervalAnnotation(vjref["TRAV2*02"], 0, 0, 0))
    cia.add(IntervalAnnotation(vjref["TRAV3*01"], 0, 0, 0))
    with pytest.raises(ValueError):
        cia.add(IntervalAnnotation(vjref["TRAJ3*01"], 0, 0, 0))
    with pytest.raises(ValueError):
        cia.add(IntervalAnnotation(trbv2, 0, 0, 0))

    # Test adding empty annotations
    cia1 = ConsensusIntervalAnnotation()
    cia2 = ConsensusIntervalAnnotation()
    cia2.add(cia1)
    assert cia2.allele is None
    assert cia2.start == 0
    assert cia2.end == 0
    assert cia2.score == 0

    # Test equality
    cia1 = ConsensusIntervalAnnotation()
    cia2 = ConsensusIntervalAnnotation()
    assert cia1 == cia2
    cia1.add(IA1)
    cia1.add(IA2)
    assert cia1 != cia2
    cia2.add(IA2)
    cia2.add(IA1)
    assert cia1 == cia2

def test_AnnotatedClone():
    from rtcr.allele import Allele
    from rtcr.clone import AnnotatedClone, CloneSet
    from rtcr.seq import IntervalAnnotation
    import cPickle as pickle
    
    gen_allele = lambda gid, seq:\
        Allele(species = "", name = gid, functionality = "F",
            refpos = 0, acc_nr = str(hash(gid)), seq = seq)

    def gen_clone(v_allele, j_allele, count = 1, oseq = "",
            score = 20):
        vseq = "" if v_allele is None else v_allele.seq
        jseq = "" if j_allele is None else j_allele.seq
        vs = 0
        ve = len(vseq)
        js = ve + len(oseq)
        je = js + len(jseq)
        seq = vseq + oseq + jseq
        v_annot = IntervalAnnotation(v_allele, vs, ve, score)
        j_annot = IntervalAnnotation(j_allele, js, je, score)
        return AnnotatedClone(v_annot, None, j_annot, None,
                (seq, [40]*len(seq), count))
    v1a1 = gen_allele("TRBV1*01", "A"*10)
    v1a2 = gen_allele("TRBV1*02", "A"*10)
    v2a1 = gen_allele("TRBV2*01", "A"*10)
    d1a1 = gen_allele("TRBD1*01", "AGCAG")
    d2a1 = gen_allele("TRBD2*01", "CCGCC")
    j1a1 = gen_allele("TRBJ1*01", "C"*10)
    j1a2 = gen_allele("TRBJ1*02", "C"*10)
    j2a1 = gen_allele("TRBJ2*01", "C"*10)
    c1a1 = gen_allele("TRBC1*01", "G"*10)
    c2a1 = gen_allele("TRBC2*01", "GGGGGCGGGG")

    c1 = gen_clone(v1a1, j1a1)
    assert c1.seq == v1a1.seq + j1a1.seq
    assert c1.v.allele == v1a1
    assert c1.j.allele == j1a1
    assert c1.count == 1

    # Adding different allele
    c2 = c1.add(gen_clone(v1a2, j1a1, 2))
    assert c2.v.allele == v1a2

    # Adding a different segment
    c2 = c2.add(gen_clone(v2a1, j1a1, 4))
    assert c2.v.allele == v2a1

    # Test adding wrong segment type
    with pytest.raises(ValueError):
        gen_clone(j1a1, v1a1)

    # Test missing V or missing J region
    with pytest.raises(ValueError):
        gen_clone(None, j1a1)
    with pytest.raises(ValueError):
        gen_clone(v1a1, None)

    # Test AnnotatedClones are distinguished based on sequence
    c1 = gen_clone(v1a1, j1a1)
    c2 = gen_clone(v2a1, j1a1)
    assert c1.v != c2.v
    assert c1.seq == c2.seq
    assert c1 == c2
    c3 = gen_clone(v1a1, j1a1, 1, "ATCAGA")
    assert c3.v == c1.v
    assert c3.d == c1.d
    # c3 has extra nucleotides, so j is in a different position
    assert c3.j != c1.j
    # however, the allele should still be the same
    assert c3.j.allele == c1.j.allele
    assert c3.c == c1.c
    assert c3.seq != c1.seq
    assert c3 != c1

    # Test adding to a CloneSet
    cs = CloneSet()
    c1 = gen_clone(v1a1, j1a1)
    c2 = gen_clone(v1a2, j1a1)
    cs.add(c1)
    assert c2 in cs # AnnotatedClone is identified by sequence alone
    with pytest.raises(KeyError):
        cs.add(c2)
    assert len(cs) == 1
    assert cs.sequence_count == 1
    cs.add(c2, merge = True) # actually add the clone
    assert len(cs) == 1
    assert cs.sequence_count == 2
    assert cs.mutation_count == 0
    c = next(iter(cs))
    assert c1.orig_count == 1
    assert c2.orig_count == 1
    assert c.orig_count == c1.orig_count + c2.orig_count 
    assert c.count == 2
    assert c.mutation_count == 0

    # Test frequency wins for allele annotation
    cs.add(gen_clone(v2a1, j1a1, 3), merge = True)
    assert len(cs) == 1
    assert cs.sequence_count == 5
    assert cs.mutation_count == 0
    c = next(iter(cs))
    assert c.v.allele == v2a1
    assert c.d.allele is None
    assert c.j.allele == j1a1
    assert c.c.allele is None
    assert c.count == 5

    cs.add(gen_clone(v2a1, j1a1, 1, "ATT"), merge = True)
    assert len(cs) == 2
    assert cs.mutation_count == 0

    # TODO: test quality wins for allele annotation (if tied on frequency)
    # TODO: test frequency wins for interval annotation
    # TODO: test quality wins for interval annotation (if tied on frequency)

    # Test pickling
    for clone in [c1, c2, c3]:
        s = pickle.dumps(clone)
        u = pickle.loads(s)
        assert not u is clone
        assert u.count == clone.count
        assert u.base_count == clone.base_count
        assert u.mutation_count == clone.mutation_count
        assert u.orig_count == clone.orig_count
        assert u.qual == clone.qual
        assert u.seq == clone.seq
        assert u.v == clone.v
        assert u.d == clone.d
        assert u.j == clone.j
        assert u.c == clone.c

def test_AnnotatedCloneDistinctAllele():
    from rtcr.allele import Allele
    from rtcr.clone import AnnotatedCloneDistinctAllele, CloneSet
    from rtcr.seq import IntervalAnnotation
    import cPickle as pickle
    
    gen_allele = lambda gid, seq:\
        Allele(species = "", name = gid, functionality = "F",
            refpos = 0, acc_nr = str(hash(gid)), seq = seq)

    def gen_clone(v_allele, j_allele, count = 1, oseq = "",
            score = 20):
        vseq = "" if v_allele is None else v_allele.seq
        jseq = "" if j_allele is None else j_allele.seq
        vs = 0
        ve = len(vseq)
        js = ve + len(oseq)
        je = js + len(jseq)
        seq = vseq + oseq + jseq
        v_annot = IntervalAnnotation(v_allele, vs, ve, score)
        j_annot = IntervalAnnotation(j_allele, js, je, score)
        return AnnotatedCloneDistinctAllele(v_annot, None, j_annot, None,
                (seq, [40]*len(seq), count))
    v1a1 = gen_allele("TRBV1*01", "A"*10)
    v1a2 = gen_allele("TRBV1*02", "A"*10)
    v2a1 = gen_allele("TRBV2*01", "A"*10)
    d1a1 = gen_allele("TRBD1*01", "AGCAG")
    d2a1 = gen_allele("TRBD2*01", "CCGCC")
    j1a1 = gen_allele("TRBJ1*01", "C"*10)
    j1a2 = gen_allele("TRBJ1*02", "C"*10)
    j2a1 = gen_allele("TRBJ2*01", "C"*10)
    c1a1 = gen_allele("TRBC1*01", "G"*10)
    c2a1 = gen_allele("TRBC2*01", "GGGGGCGGGG")

    # Test generating the clone
    c1 = gen_clone(v1a1, j1a1)
    assert c1.seq == v1a1.seq + j1a1.seq
    assert c1.v.allele == v1a1
    assert c1.j.allele == j1a1
    assert c1.count == 1

    # Test clones are not equal when annotated with different allele or segment
    c2 = gen_clone(v1a1, j1a1)
    c3 = gen_clone(v1a2, j1a1) # different v allele
    c4 = gen_clone(v2a1, j1a1) # different v segment
    c5 = gen_clone(v1a1, j1a2) # different j allele
    c6 = gen_clone(v1a1, j2a1) # different j segment
    assert c2 == c1
    assert c3 != c2
    assert c4 != c3
    assert c4 != c1
    assert c5 != c1
    assert c6 != c1

def helperfunc_test_get_error_stats(load_c):
    if load_c:
        from cseq import get_error_stats
    else:
        from rtcr.seq import get_error_stats

    from rtcr.align import SAMRecord

    gen_sam = lambda seq, qual, cigar, ras:\
            SAMRecord(QNAME = "", FLAG = 0, RNAME = "", POS = ras + 1,
                    MAPQ = 0, CIGAR = cigar, RNEXT = '*', PNEXT = 0,
            TLEN = 0, SEQ = seq, QUAL = "".join([chr(q + 33) for q in qual]))
    q_mm = [0] * 42
    q_n = [0] * 42

    # Test simple base counting without any sequence errors
    sam1 = gen_sam("ACTTG", [40]*5, "5M", 0)
    n, mm, ins, dels, r_roi_as, r_roi_ae = get_error_stats(sam1, "ACTTG",
            q_mm, q_n)
    assert n == 5
    assert mm == 0
    assert ins == 0
    assert dels == 0
    assert q_mm[40] == 0
    assert q_n[40] == 5
    assert r_roi_as == 0
    assert r_roi_ae == 5

    # Test q_mm and q_n are updated
    n, mm, ins, dels, r_roi_as, r_roi_ae = get_error_stats(sam1, "ATTTG",
            q_mm, q_n)
    assert n == 5
    assert mm == 1
    assert ins == 0
    assert dels == 0
    for i in xrange(42):
        if i != 40:
            assert q_mm[i] == 0
            assert q_n[i] == 0
        else:
            assert q_mm[i] == 1
            q_n[i] == 10
    assert r_roi_as == 0
    assert r_roi_ae == 5

    # Test base counting with varying quality scores, indels, and mismatches
    # ACTTTAAGGACCCGCTTG---T (reference seq)
    #   SSS=X=DD===X====III= (expanded cigar showing matches and mismatches)
    #   GAAATG--CCCACTTGGCCT (query)
    rseq = "ACTTTAAGGACCCGCTTGT"
    sam2 = gen_sam("GAAATGCCCACTTGGCCT",
            # phred in softclip, these should not be counted
            [41] * 3 + \
            [40] * 6 + \
            # second mismatch gets different phred score
            [28] + \
            [40] * 8, "3S3M2D8M3I1M", 5)
    n, mm, ins, dels, r_roi_as, r_roi_ae = get_error_stats(sam2, rseq,
            q_mm, q_n)
    assert n == 15
    assert mm == 2
    assert ins == 3
    assert dels == 2
    assert q_mm[28] == 1
    assert q_mm[40] == 2
    assert q_mm[41] == 0
    assert q_n[41] == 0 # (soft-clip test)
    assert q_n[40] == 24
    assert r_roi_as == 5
    assert r_roi_ae == 19

    # Test if "=" and "X" characters are ignored in favor of checking against
    # reference sequence.
    q_mm2 = [0]*42
    q_n2 = [0]*42
    sam3 = gen_sam("AAA", [40]*3, "2X1=", 0)
    n, mm, ins, dels, r_roi_as, r_roi_ae = get_error_stats(sam3, "AAA", q_mm2,
            q_n2)
    assert n == 3
    assert mm == 0
    assert ins == 0
    assert dels == 0
    assert q_mm2[40] == 0
    assert q_n2[40] == 3
    assert r_roi_as == 0
    assert r_roi_ae == 3
    n, mm, ins, dels, r_roi_as, r_roi_ae = get_error_stats(sam3, "TTT", q_mm2,
            q_n2)
    assert n == 3
    assert mm == 3
    assert ins == 0
    assert dels == 0
    assert q_mm2[40] == 3
    assert q_n2[40] == 6
    assert r_roi_as == 0
    assert r_roi_ae == 3

    # Test getting error stats for a region of interest only
    q_mm3 = [0]*42
    q_n3 = [0]*42
    sam4 = gen_sam("A"*5 + "C"*5, [40]*5 + [30]*5, "10M", 0)
    n, mm, ins, dels, r_roi_as, r_roi_ae = get_error_stats(sam4, "AXAAACCTTC",
            q_mm3, q_n3, 4, 8)
    assert n == 4
    assert mm == 1
    assert q_mm3[30] == 1
    assert q_mm3[40] == 0
    assert q_n3[40] == 1
    assert q_n3[30] == 3
    assert ins == 0
    assert dels == 0
    assert r_roi_as == 0
    assert r_roi_ae == 4

    # Test getting error stats for a region of interest containing indels
    q_mm4 = [0] * 42
    q_n4 = [0] * 42
    sam5 = gen_sam("AAACCCGGGAAGGG", [10]*3 + [20]*3 + [30]*3 + [40]*5,
            "2S3M5D4M2I3M",0)
    n, mm, ins, dels, r_roi_as, r_roi_ae = get_error_stats(sam5,
            #AAACC-----CGGGAAGGG
            #   *        **II * 
              "ATCGGGGGCGCCGTG", q_mm4, q_n4, 2, 10)
    assert n == 3
    assert mm == 0
    assert ins == 0
    assert dels == 5
    for i in xrange(42):
        assert q_mm4[i] == 0
    assert q_n4[20] == 2
    assert q_n4[30] == 1
    assert r_roi_as == 0
    assert r_roi_ae == 8
    q_mm4 = [0] * 42
    q_n4 = [0] * 42
    n, mm, ins, dels, r_roi_as, r_roi_ae = get_error_stats(sam5,
            "ATCGGGGGCGCCGTG", q_mm4, q_n4, 8, 12)
    assert n == 4
    assert mm == 2
    assert ins == 0
    assert dels == 0
    assert q_mm4[20] == 0 
    assert q_mm4[30] == 2
    assert q_mm4[40] == 0 
    assert q_n4[20] == 1
    assert q_n4[30] == 3
    assert q_n4[40] == 0
    assert r_roi_as == 0
    assert r_roi_ae == 4

    # Test where alignment starts and ends within the ROI
    Q_mm = [0] * 42
    Q_n = [0] * 42
    #   _________________________  (ROI)
    # CGTGGGCGG-TGGA--GGCTGGTCGCGA (reference)
    #       SS   *
    #       GGGATCGATTGG-TGG       (query)
    sam6 = gen_sam("GGGATCGATTGGTGG", [40]*5 + [20] + [40]*10,
            "2S1M1I4M2I2M1D3M", 8)
    n, mm, ins, dels, r_roi_as, r_roi_ae = get_error_stats(sam6,
            "CGTGGGCGGTGGAGGCTGGTCGCGA", Q_mm, Q_n, 2, 23)
    assert n == 13
    assert mm == 1
    assert ins == 3
    assert dels == 1
    assert Q_mm[20] == 1
    assert Q_n[20] == 1
    assert Q_mm[40] == 0
    assert Q_n[40] == 12
    assert r_roi_as == 6
    assert r_roi_ae == 17

def test_get_error_stats_c():
    helperfunc_test_get_error_stats(load_c = True)

def test_get_error_stats_python():
    helperfunc_test_get_error_stats(load_c = False)

def test_build_clone(vjref):
    from rtcr.align import SAMRecord
    from rtcr.clone import AnnotatedClone, build_clone

    random.seed(1)

    # Create sequence error free clone
    preseq = randseq(21)
    vid = "TRAV5*03"
    jid = "TRAJ2*02"
    vseq = vjref[vid].seq
    oseq = randseq(10)
    jseq = vjref[jid].seq
    seq = preseq + vseq + oseq + jseq
    qual = 'I' * len(seq)
    seq2 = seq[len(preseq) + len(vseq):]
    qual2 = qual[len(preseq)+ len(vseq):]
    nt_jnc_v = vseq[vjref[vid].refpos:]
    nt_jnc_j = jseq[:vjref[jid].refpos]
    nt_jnc = nt_jnc_v + oseq + nt_jnc_j
#    print seq
#    print " "*len(preseq) + "v"*len(vseq) + "o"*len(oseq) + "j"*len(jseq)
#    print " "*(len(preseq) + len(vseq) - len(nt_jnc_v)) + nt_jnc
    v_rec = SAMRecord(QNAME = "", FLAG = 0, RNAME = vid, POS = 1, MAPQ = 9,
            CIGAR = '%sS%sM'%(len(preseq), len(vseq)), RNEXT = '*', PNEXT = 0,
            TLEN = 0, SEQ = seq, QUAL='I' * len(seq))
    j_rec = SAMRecord(QNAME = "", FLAG = 0, RNAME = jid, POS =1, MAPQ =17,
            CIGAR='%sS%sM'%(len(oseq), len(jseq)), RNEXT = '*', PNEXT = 0,
            TLEN = 0, SEQ = seq2,
            QUAL = 'I' * len(seq2))

    
    c1 = build_clone(vjref, v_rec, j_rec, False)
    assert c1.v.allele.name == vid
    assert c1.j.allele.name == jid
    assert c1.seq == nt_jnc

    c2 = build_clone(vjref, v_rec, j_rec, False, "Clone")
#    print c2
#    assert False

def test_get_vj_alignments(vjref):
    from rtcr.util import TemporaryDirectory
    from rtcr.fileio import FastqFormat
    from rtcr.seq import QSequence
    from rtcr.align import get_vj_alignments
    from itertools import cycle, izip, count
    import os

    random.seed(1)

    v_alleles = vjref.get_alleles(region = "V")
    j_alleles = vjref.get_alleles(region = "J")

    # Prepare test dataset
    qseqs = []
    for i, v_allele, j_allele in izip(
            count(),
            cycle(v_alleles),
            cycle(j_alleles)):
        if i >= 1000:
            break
        # in the test reference the V and J alleles are guaranteed to not end
        # on a "T" nucleotide. Hence, here we add the "T" to properly mark the
        # start of the region between V and J.
        oseq = "T" + randseq(10) + "T"
        pre_seq = randseq(random.randint(0,20)) + "T"
        full_seq = pre_seq + v_allele.seq + oseq + j_allele.seq
        qseq = QSequence("%s|%s|%s|%s|"%(v_allele.name, oseq, j_allele.name,
            len(pre_seq)),
                full_seq, [40]*len(full_seq))
        if i % 2 == 1:
            qseq = qseq.reverse_complement(qseq.name + "rc")
        qseqs.append(qseq)
    assert len(qseqs) == 1000

    with TemporaryDirectory() as tmpdir:
        fq_fn = os.path.join(tmpdir, "ref.fq")
        FastqFormat.records_out(open(fq_fn, 'w'),
                qseqs)
        cmd_build_index = "bowtie2-build"
        args_build_index = r"%(ref_fn)s %(index_fn)s"
        cmd_align = "bowtie2"
        args_align_base = r"--local -x %(index_fn)s - " + \
                "--phred%(phred_encoding)s --threads %(n_threads)s"
        args_align_v = r"-D 20 -R 3 -N 0 -i S,1,0.50 -L 8 " + args_align_base
        args_align_j = r"-D 20 -R 3 -N 0 -i S,1,0.50 -L 6 " + args_align_base
        results = {v_rec.QNAME:(v_rec, j_rec) \
                for v_rec, j_rec in get_vj_alignments(vjref, open(fq_fn, 'r'),
                cmd_build_index, 
                args_build_index, cmd_align, args_align_v, args_align_j,
                phred_encoding = 33, n_threads = 7)}
        assert len(results) == 1000
        for qseq in qseqs:
            assert qseq.name in results

            v_rec, j_rec = results[qseq.name]
            assert v_rec.QNAME == j_rec.QNAME
            vid, oseq, jid, len_pre_seq, rc = v_rec.QNAME.split("|")
            len_pre_seq = int(len_pre_seq)
            assert v_rec.RNAME == vid
            assert v_rec.CIGAR == "%sS50M62S"%len_pre_seq
            assert j_rec.RNAME == jid
            assert j_rec.CIGAR == "12S50M"

#def test_qmerge():
#    from rtcr.fileio import zopen, LegacyCloneSetFormat 
#    from rtcr.clone import Clone, CloneSet
#    from qmerge import run_qmerge_on_bin
#
#    gen_clone = lambda seq,qual,count:Clone((seq, qual, count))
#    c0 = gen_clone("AAAA", [40, 40, 40, 40], 10)
#    c1 = gen_clone("AATA", [40, 40, 1, 40], 1)
#    c2 = gen_clone("AATA", [40, 40, 1, 40], 1)
#
#    run_qmerge_on_bin(CloneSet((c0, c1)), mismatch_rate = .25, confidence = .5,
#            max_Q = 30)
#    # TODO: find smallish dataset to run qmerge on
#    cloneset = LegacyCloneSetFormat.records_in(
#            zopen("/home/bram/tmp/hs25_10x_cnd_raw_len45.tsv", 'r'))
#    run_qmerge_on_bin(cloneset, mismatch_rate = 0.0106028457443,
#            confidence = .5, max_Q = 30)
#
#def test_imerge():
#    from itertools import chain
#    from rtcr.clone import Clone, CloneSet
#    from imerge import run_imerge
#
#    c0 = Clone(("A"*45, [40]*45, 65))
#    c1 = Clone(("A"*44 + "T", [40]*45, 29))
#    c2 = Clone(("A"*43 + "TT", [40]*45, 5))
#    c3 = Clone(("A"*42 + "TTT", [40]*45, 1))
#
#    cloneset = CloneSet()
#    cloneset.add(c0)
#    cloneset.add(c1)
#    cloneset.add(c2)
#    cloneset.add(c3)
#
#    cloneset2 = CloneSet(cloneset)
#    assert cloneset2.count == 4
#    assert cloneset2.sequence_count == 100
#    assert cloneset2.base_count == 4500
#    assert cloneset2.mutation_count == 0
#
#    cloneset = run_imerge(cloneset, mismatch_rate = .01, confidence = .99)
#    assert cloneset.count == 1
#    assert cloneset.sequence_count == 100
#    assert cloneset.base_count == 4500
#    assert cloneset.mutation_count == 42 
#
#    cloneset3 = CloneSet()
#    cloneset3.add(Clone(("AAA", [40]*3, 10)))
#
#    cloneset4 = CloneSet(chain.from_iterable([cloneset2, cloneset3]))
#    assert cloneset4.count == 5
#
#def test_imerge_equivalence():
#    #Test re-implementation of RTCR's imerge algorithm, to see if it is
#    # equivalent to the old (version 0.2.2).
#    from rtcr.clone import Clone, CloneSet, SearchableCloneSet
#    from rtcr.fileio import zopen, LegacyCloneSetFormat 
#    from imerge import run_imerge 
#
#    # Load example dataset
#    cloneset = LegacyCloneSetFormat.records_in(
#            zopen("/home/bram/tmp/hs25_10x_cnd_raw_qmerge_len45.tsv", 'r'))
#
#    #print "#clones: %s"%len(cloneset)
#    #scls = SearchableCloneSet(cloneset)
#    
#    #correct_neighbors = [\
#    #        (0, "TGTGCCAGCAGTTACTCCCGGGGCGGGGGCCCCACGCAGTATTTT", 147),\
#    #        (2, "TGTGCCAGCAGTAACTCCCGGGGCGGGGGCCCCACGCAGTGTTTT", 1)]
#    #assert set([(hd, clone.seq, clone.count) for hd, clone in \
#    #        scls.neighbors("TGTGCCAGCAGTTACTCCCGGGGCGGGGGCCCCACGCAGTATTTT",
#    #    maxhd=2)]) == set(correct_neighbors)
#    run_imerge(cloneset, mismatch_rate = 0.0106028457443,
#            confidence = .99)
