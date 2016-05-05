# Description: provides datastructures to hold (DNA) sequences.

from string import maketrans
from collections import deque, namedtuple
from itertools import izip, islice, ifilter
import re
from operator import itemgetter

DNA_COMPLEMENT = maketrans("ATCG","TAGC")

class Sequence(object):
    """Immutable sequence."""
    def __init__(self, name, seq):
        assert(isinstance(name, basestring))
        assert(isinstance(seq, basestring))
        self._name = name
        self._seq = seq
        self._hash = hash((self._name, self._seq))

    def reverse_complement(self, name = ""):
        """Return reverse complement as Sequence object."""
        return Sequence(name = name, \
                seq = self._seq.translate(DNA_COMPLEMENT)[::-1])

    @property
    def name(self):
        return self._name

    @property
    def seq(self):
        return self._seq

    def __add__(self, other):
        return Sequence(self._seq + other._seq)

    def __len__(self):
        return len(self._seq)

    def __eq__(self, other):
        if isinstance(other, Sequence):
            return self._name == other._name and self._seq == other._seq
        else:
            return NotImplementedError()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return self._hash

    def __repr__(self):
        return "%s(name: \"%s\", seq: \"%s\")%"%(self.__class__.__name__,
                self._name, self._seq)

    def __str__(self):
        return self._seq

class QSequence(Sequence):
    """Immutable Sequence with for each base a (Phred) quality score."""

    def __init__(self, name, seq, qual, force_N_to_zero = False):
        super(QSequence, self).__init__(name, seq)
        assert(len(seq) == len(qual))

        if force_N_to_zero:
            for i in xrange(len(seq)):
                if seq[i] == 'N':
                    qual[i] = 0

        self._qual = tuple(qual)
        self._hash = hash((self._name, self._seq, self._qual))

    def reverse_complement(self, name = ""):
        return QSequence(name = name,
                seq = super(QSequence, self).reverse_complement(name).seq,
                qual = self._qual[::-1])

    @property
    def qual(self):
        return self._qual

    def __add__(self, other):
        return QSequence("(%s,%s)"%(self._name, other._name),
                self._seq + other._seq, self._qual + other._qual)

    def __eq__(self, other):
        if isinstance(other, QSequence):
            return self._name == other._name and self._seq == other._seq and \
                    self._qual == other._qual
        else:
            raise NotImplementedError()

    def __repr__(self):
        return "%s(name: \"%s\", seq: \"%s\", qual: \"%s\")"%(
                self.__class__.__name__, self._name, self._seq, self._qual)

    def __str__(self):
        return "%s:%s"%(self._seq, self._qual)

class PfmElement(object):
    """Element of a position frequency matrix, containing base frequencies and
    their maximum observed quality score.
    """
    def __init__(self):
        self._count = 0
        self._base = "N"
        self._q = 0

        self._max_resolved_freq = 0
        # the lambda in the defaultdict makes objects unpicklable
        #self._basefreqs = defaultdict(lambda:[0,0])
        self._basefreqs = {}

    @property
    def base(self):
        return self._base

    @property
    def q(self):
        return self._q

    @property
    def count(self):
        """Return total number of bases (including unresolved bases)."""
        return self._count

    @property
    def mutation_count(self):
        """Return number of mutations required to turn all bases into the most
        common resolved base.
        """
        if self._max_resolved_freq > 0:
            return self.count - self._max_resolved_freq
        else:
            return 0

    def add(self, other):
        """Add other PfmElement."""
        for base, (freq, q) in other._basefreqs.iteritems():
            self.inc_freq(base, q, freq)

    def inc_freq(self, base, q, count):
        """Increase frequency of a particular base."""
        assert count > 0
        assert base != "N" or q == 0
        if not base in self._basefreqs:
            self._basefreqs[base] = [0, 0]
        self._basefreqs[base][0] += count
        self._count += count
        self._basefreqs[base][1] = max(self._basefreqs[base][1], q)

        if self._basefreqs[base][0] < self._max_resolved_freq or \
                base == "N":
            return

        if self._basefreqs[base][0] > self._max_resolved_freq:
            self._max_resolved_freq = self._basefreqs[base][0]
            self._base = base
            self._q = self._basefreqs[base][1]
            return

        self._max_resolved_freq = self._basefreqs[base][0]

        # Break tie by selecting base with highest quality score.
        tied = sorted(\
                ifilter(lambda (base, (freq, q)): base != "N" and \
                freq == self._max_resolved_freq, self._basefreqs.iteritems()),
                        key = lambda (base, (freq, q)): -q)

        # Structure for elements of "tied": (base, [freq, q]).
        # Example for tied: [('C', [10, 30]), ('T', [10, 30])]
        if tied[0][1][1] == tied[1][1][1]: # q scores are equal
            self._base = "N"
            self._q = 0
        else:
            self._base = tied[0][0]
            self._q = tied[0][1][1]

    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__,
                ";".join(["%s:%s,Q%s"%(base, freq, q)\
                for base, (freq, q) in self._basefreqs.iteritems()]))

    def __str__(self):
        return ";".join(["[%s%sQ%s]"%(base, freq, q) if base == self.base \
                else "%s%sQ%s"%(base, freq, q) \
                for base, (freq, q) in self._basefreqs.iteritems()])

class ConsensusQSequence(object):
    """Contains consensus sequence and quality scores.

    The consensus at every position is the most common base with the highest
    observed quality score for that base. In case of ties on frequency, an "N"
    is inserted with a quality score of 0.
    """
    def __init__(self, other = None):
        if other is None:
            seq = ""
            qual = []
            count = 0
        elif isinstance(other, ConsensusQSequence):
            self._pfm = deque()
            self._seq = ""
            self._qual = []
            self._count = 0
            self._base_count = 0
            self._mutation_count = 0
            self.add(other);
            return
        elif hasattr(other, "seq"):
            seq = other.seq
            qual = list(other.qual)
            if hasattr(other, "count"):
                count = other.count
            else:
                count = 1
        elif isinstance(other, tuple) and len(other) == 3 and \
                isinstance(other[0], tuple):
            # Filling pfm based on a tuple, assuming alphabet is "ACGTN"
            freq, qual, count = other
            if len(freq) != len(qual):
                raise ValueError("freq and qual in tuple (first two items) \
are not the same length")
            if len(freq) % 5 != 0:
                raise ValueError("length of freq (first item in tuple) is not \
a multiple of 5")
            if count > sum(freq):
                raise ValueError("number of bases in pfm exceeds sequence \
count.")
            pfm = deque()
            self._pfm = deque()
            self._seq = ""
            self._qual = []
            self._count = 0 
            self._base_count = 0
            self._mutation_count = 0
            for i, (f, q) in enumerate(izip(freq, qual)):
                if i % 5 == 0:
                    el = PfmElement()
                    pfm.append(el)
                if f > 0:
                    el.inc_freq("ACGTN"[i % 5], q, f)
            self.add(pfm);
            self._count = count
            return
        else:
            seq, qual = other[:2]
            if len(other) == 2:
                count = 1
            else:
                count = other[2]

        assert count >= 0
        assert len(seq) == len(qual)

        pfm = deque()
        for base, q in izip(seq, qual):
            assert(0 <= q <= 41)
            assert(base in "ACGTN")
            el = PfmElement()
            el.inc_freq(base, q, count)
            pfm.append(el)

        self._seq = seq
        self._qual = qual
        self._pfm = pfm 
        self._count = count
        self._base_count = count * len(seq) 
        self._mutation_count = 0

    @property
    def seq(self):
        return self._seq

    @property
    def qual(self):
        return self._qual

    @property
    def count(self):
        """Return the number of sequences."""
        return self._count

    @property
    def mutation_count(self):
        """Return the number of mutated bases."""
        return self._mutation_count

    @property
    def base_count(self):
        """Return the total number of bases."""
        return self._base_count

    @property
    def pfm(self):
        return self._pfm

    def _accommodate(self, offset, seq_len):
        # Example for "offset":
        # consensus = ACTG
        # other.seq = AACTG
        #
        # offset = -1 because "other.seq" has to be moved 1 position to the
        # left to align with the consensus(/pfm).

        for i in xrange(-offset):
            self._pfm.appendleft(PfmElement())
        self._seq = "N"*(-offset) + self._seq
        self._qual = [0]*(-offset) + self._qual

        if offset < 0:
            offset = 0

        n = offset + seq_len - len(self._pfm)

        for i in xrange(n):
            self._pfm.append(PfmElement())

        self._seq += "N"*n
        self._qual += [0]*n

        return offset

    def _update(self, pos, pfm_el, other):
        self._base_count -= pfm_el.count
        self._mutation_count -= pfm_el.mutation_count

        if isinstance(other, PfmElement):
            pfm_el.add(other)
        else:
            pfm_el.inc_freq(*other)

        self._base_count += pfm_el.count
        self._mutation_count += pfm_el.mutation_count

        if self._seq[pos] != pfm_el.base:
            self._seq = self._seq[:pos] + pfm_el.base + self._seq[pos+1:]

        self._qual[pos] = pfm_el.q
        self._count = max(self._count, pfm_el.count)

    def add(self, other, offset = 0):
        """Add ConsensusQSequence, position frequency matrix (pfm),
        or a QSequence(-like) object.

        :other: a pfm (should be derived from deque), or an object containing
        "pfm" as attribute. If it is not a pfm, it should be an object
        containing "seq" and "qual" attributes containing a sequence and
        corresponding quality scores respectively.
        :offset: number of bases the other sequence needs to be moved to the
        right in order to align to the consensus sequence.
        """
        if (offset > 0 and offset > len(self)) or \
            (offset < 0 and abs(offset) > len(other)):
            raise ValueError("offset so large/small that there is a gap")
        offset = self._accommodate(offset, len(other))
        pfm = None
        if isinstance(other, deque):
            pfm = other
        elif hasattr(other, "pfm"):
            pfm = other.pfm

        if hasattr(other, "count") and isinstance(other.count, int): 
            new_count = self._count + other.count
        else:
            new_count = None
        
        if not pfm is None:
            for i, (self_el, el) in enumerate(izip(
                    islice(self._pfm, offset, None), pfm)):
                self._update(i + offset, self_el, el)
        else:
            for i, (self_el, base, q) in \
                    enumerate(izip(
                    islice(self._pfm, offset, None), other.seq, other.qual)):
                self._update(i + offset, self_el, (base, q, 1))

        if not new_count is None:
            self._count = new_count

    def __add__(self, other):
        self.add(other)

    def __len__(self):
        return len(self._seq)

    def __repr__(self):
        return "%s(seq: \"%s\", qual: %s, pfm: [%s])"%(
                self.__class__.__name__,
                self._seq,
                self._qual,
                "|".join([repr(el) for el in self._pfm]))

    def __str__(self):
        return "|".join(map(str, [el for el in self._pfm]))

def merge_stats(cqs1, cqs2, offset = 0):
    if offset < 0:
        cqs1, cqs2 = cqs2, cqs1
        offset = -offset

    hd = 0
    merge_score = 0
    max_mut_q = 0

    len_cqs1 = len(cqs1.seq)
    len_cqs2 = len(cqs2.seq)
    if len_cqs1 > 0 and len_cqs2 > 0:
        for i in xrange(max(len_cqs1, len_cqs2 + offset)):
            if cqs1.seq[i] != cqs2.seq[i + offset]:
                minq = min(cqs1.qual[i], cqs2.qual[i + offset])
                merge_score += minq
                max_mut_q = max(max_mut_q, minq)
                hd += 1
    return hd, merge_score, max_mut_q

def get_offset(s1, s2, k, moff):
    """Return number of positions s1 should be moved to the right (or s2 to the
    left) so as to maximize the overlap between the two sequences. Returns
    2-tuple (n, m), where n is the offset and m is the overlap score (see
    below).

    Note, an offset of 2 means s1 should be moved 2 positions to the right.
    (e.g. s1 =   "ACTG",
          s2 = "AAACTG").

    It selects a window (of size k) in the middle of s1 and compares that
    with k-mers at the same coordinate of s2, but shifts this coordinate by
    no more than "moff" bases. The offset providing highest overlap score is
    returned. The overlap score is calculated by counting -1 for mismatches, 1
    for matches, and 0 for any (mis)matches involving "N" bases.

    :s1: DNA string
    :s2: DNA string
    :k: size of the window in s1 and s2 that is used to compare the strings.
    :moff: maximum offset of the window, i.e. number of positions to shift the
    window by to the left and to the right.
    """
    swapped = False
    # Make sure s2 is the longer of the two strings because different
    # kmers are going to be selected within s2.
    if len(s1) > len(s2):
        swapped = True
        s1,s2 = s2,s1

    # Get midpoint
    mid = len(s1) // 2

    # Make sure k-mer is not longer than the strings
    if k > min(len(s1),len(s2)):
        raise ValueError("String(s) shorter than window")

    # Get kmer from s1
    start = mid - ( k // 2 )
    end = start + k
    kmer_s1 = s1[start : end]

    # Make sure no offsets are visited that would select a kmer outside of
    # s2.
    moff_left = min(moff, start)
    moff_right = min(moff, len(s2) - end)

    # Get offset that provides most overlap
    best_match = None
    for abs_offset in xrange(max(moff_left, moff_right) + 1):
        for flip, moff in ((-1, moff_left), (1, moff_right)):
            if (flip == 1 and abs_offset == 0) or abs_offset > moff:
                continue

            offset = abs_offset * flip
            kmer_s2 = s2[start + offset : end + offset]

            # Calculate number of matches
            m = 0
            for ch1, ch2 in izip(kmer_s1, kmer_s2):
                if ch1 == 'N' or ch2 == 'N':
                    m += 0
                elif ch1 == ch2:
                    m += 1
                else:
                    m -= 1

            if best_match is None or m > best_match[1]:
                best_match = (offset, m)
            if m == k: # full match found
                return (-best_match[0], best_match[1]) if swapped else \
                        (best_match[0], best_match[1])
    return (-best_match[0], best_match[1]) if swapped else \
            (best_match[0], best_match[1])

CIGARCHARS = "MIDNSHP=X"
RE_CIGARCHARS = "[" + CIGARCHARS + "]"

def cigar_get_numbers(cigar):
    return map(int, re.split(RE_CIGARCHARS, cigar)[:-1])

def cigar_get_chars(cigar):
    return re.findall(RE_CIGARCHARS, cigar)

def cigar_intervals(cigar, ras):
    """Find start and end of the alignment intervals on the query and
    reference sequences, using the CIGAR string. Returns 4-tuple
    (qas, qae, ras, rae), where qas and qae are the 0-based alignment
    start- and end-positions on the query respectively (SEQ field of a SAM
    record), followed by a similar pair for the reference.

    :ras: 0-based(!) lefmost mapping position on reference. Note, the 'POS'
    field of a SAM record is 1-based, therefore subtract 1 from that field
    when passing it to this method.
    """
    cigar_numbers = cigar_get_numbers(cigar)
    cigar_chars = cigar_get_chars(cigar)

    rae = ras
    qas = 0
    qae = qas

    for i in xrange(len(cigar_chars)):
        n = cigar_numbers[i]
        ch = cigar_chars[i]
        if ch in "MX=":
            rae += n
            qae += n
        elif ch == "I": # insertion to the reference
            qae += n
        elif ch == "D": # deletion from the reference
            rae += n
        elif ch == "N": # skipped region from the reference
            rae += n
        elif ch == "S": # Soft clip on the query 
            if i == 0 or (i == 1 and cigar_chars[0]=="H"):
                qas += n
                qae += n
        # No need to handle "P" as it is for properly displaying
        # multi-alignments

    return qas, qae, ras, rae

def cigar_rpos2qpos(cigar, ras, target_rpos):
    """Find which position on the query corresponds to the given position
    on the reference.

    :cigar: cigar string from a SAM record
    :ras: 0-based alignment start position on the reference
    :target_rpos: target position on the reference
    """
    if target_rpos < ras:
        return None
    cigar_numbers = cigar_get_numbers(cigar)
    cigar_chars = cigar_get_chars(cigar)

    #lm = None # last position on query that matched to the reference
    rpos = ras
    qas = 0
    qpos = qas

    for i in xrange(len(cigar_chars)):
        n = cigar_numbers[i]
        ch = cigar_chars[i]
        if ch in "MX=":
            rpos += n
            qpos += n
            #lm = qpos - 1
        elif ch == "I": # insertion to the reference
            qpos += n
        elif ch == "D": # deletion from the reference
            rpos += n
        elif ch == "N": # skipped region from the reference
            rpos += n
        elif ch == "S": # Soft clip on the query 
            qpos += n
        # No need to handle "P" as it is for properly displaying
        # multi-alignments
        if target_rpos < rpos:
            if ch in "MX=":
                return qpos - (rpos - target_rpos)
            else:
                #return lm
                return None
    return None

#def map_roi(cigar, qseq, rseq, ras, r_roi_start, r_roi_end, Q_mm, Q_n):
#    """Finds where alignment of the region of interest on the reference starts
#    and ends on the query. Returns 7-tuple (roi_qas roi_qae, roi_ras,
#    roi_rae, mm, ins, dels), where roi_qas is the alignment start within the
#    region of interest as it is mapped onto the query.
#
#    :cigar: CIGAR string as defined by SAM format specification. Should contain
#    alignment of query (qseq) against reference sequence (rseq)
#    :qseq: query sequence
#    :rseq: reference sequence
#    :ras: 0-based(!) lefmost mapping position on reference. Note, the 'POS'
#    field of a SAM record is 1-based, therefore subtract 1 from that field
#    when passing it to this method.
#    :r_roi_start: 0-based start of region of interest on the reference sequence
#    :r_roi_end: end of region of intererest on the reference sequence
#    :Q_mm: list, containing for every Phred score the number of mismatches
#    (found so far).
#    :Q_n: list, containing for every Phred score the number of bases
#    (encountered so far).
#    """
#    cigar_numbers = CIGAR.get_numbers(cigar)
#    cigar_chars = CIGAR.get_chars(cigar)
#
#    rae = ras
#    qas = 0
#    qae = qas
#
#    for i in xrange(len(cigar_chars)):
#        n = cigar_numbers[i]
#        ch = cigar_chars[i]
#        if ch in "MX=":
#            rae += n
#            qae += n
#        elif ch == "I": # insertion to the reference
#            qae += n
#        elif ch == "D": # deletion from the reference
#            rae += n
#        elif ch == "N": # skipped region from the reference
#            rae += n
#        elif ch == "S": # Soft clip on the query 
#            if i == 0 or (i == 1 and cigar_chars[0]=="H"):
#                qas += n
#                qae += n
#        # No need to handle "P" as it is for properly displaying
#        # multi-alignments
#
#    return qas, qae, ras, rae

def get_error_stats(aln, rseq, Q_mm, Q_n, r_roi_start = 0, r_roi_end = 0):
    """Determines number of mismatches, insertions, and deletions in an
    alignment. Returns 6-tuple (qlen, mm, ins, dels, r_roi_as, r_roi_ae) where:
        - n, number of query bases inside the alignment (i.e. not counting
        soft/hard-clipped bases at the flanks)
        - mm, the number of mismatches
        - ins, the number of insertions
        - dels, the number of deletions
        - r_roi_as, alignment start within the region of interest
        - r_roi_ae, alignment end within the region of interest
    The Q_mm parameter is updated by adding the number of mismatches for each
    Phred score, and the Q_n parameter is updated by adding the number of bases
    that had each particular Phred score.

    If the alignment spanned the whole region of interest on the reference,
    then r_roi_as = 0 and r_roi_ae = r_roi_end - r_roi_start (i.e. the length
    of the region of interest). All the statistics (n, mm, ins, dels) including
    the updates to Q_mm and Q_n are confined to the part of the alignment
    inside the region of interest given by r_roi_start and r_roi_end.

    :aln: SAMRecord-like object containing an alignment
    :rseq: reference sequence belonging to the alignment
    :Q_mm: list-like object, indexed by Phred score, holding the number of
    mismatches (found so far) for every Phred score.
    :Q_n: list-like object, indexed by Phred score, holding the number of bases
    (found so far) for every Phred score.
    :r_roi_start: 0-based start of region of interest on the reference.
    :r_roi_end: end of region of interest on the reference. 0 means until the
    end of the reference.
    """
    if r_roi_end == 0:
        r_roi_end = len(rseq)
    if r_roi_start < 0:
        raise ValueError("roi starting before position 0")
    if r_roi_start >= r_roi_end:
        raise ValueError("roi starts after its end")

    qlen = 0
    mm = 0
    ins = 0
    dels = 0

    cigar_numbers = cigar_get_numbers(aln.CIGAR)
    cigar_chars = cigar_get_chars(aln.CIGAR)

    rpos = aln.POS - 1 # start of alignment on reference
    r_roi_as = rpos - r_roi_start if rpos > r_roi_start else 0
    r_roi_ae = r_roi_as
    qpos = 0 # position on query
    qseq = aln.SEQ
    qqualstr = aln.QUAL # ascii encoded (Phred + 33) quality string

    for i in xrange(len(cigar_chars)):
        n = cigar_numbers[i]
        ch = cigar_chars[i]
        if rpos >= r_roi_end:
            break
        if ch in "MX=":
            for j in xrange(n):
                if rpos < r_roi_start:
                    qpos += 1
                    rpos += 1
                    continue
                if rpos >= r_roi_end:
                    break
                qlen += 1
                q = ord(qqualstr[qpos]) - 33
                Q_n[q] += 1
                if rseq[rpos] != qseq[qpos]:
                    mm += 1
                    Q_mm[q] += 1
                qpos += 1
                rpos += 1
                r_roi_ae += 1
        elif ch == "I": # insertion to the reference
            qpos += n
            if rpos > r_roi_start:
                ins += n
                qlen += n
                for j in xrange(n):
                    q = ord(qqualstr[qpos]) - 33
                    Q_n[q] += 1
        elif ch == "D": # deletion from the reference
            for j in xrange(n):
                if rpos < r_roi_start:
                    rpos += 1
                    continue
                if rpos > r_roi_start and rpos < r_roi_end:
                    dels += 1
                    r_roi_ae += 1
                rpos += 1
        elif ch == "N": # skipped region from the reference
            rpos += n
        elif ch == "S": # Soft clip on the query 
            qpos += n

    return (qlen, mm, ins, dels, r_roi_as, r_roi_ae)

IntervalAnnotation = namedtuple("IntervalAnnotation", \
        ["allele", "start", "end", "score"])

#class IntervalAnnotation(object):
#    """Holds annotation for an interval on a sequence"""
#    def __init__(self, allele = None, start = 0, end = 0, score = 0):
#        if not start >= 0:
#            raise ValueError("start should be >= 0")
#        if not end >= 0:
#            raise ValueError("end should be >= 0")
#        if not end >= start:
#            raise ValueError("end should be >= 0")
#        if not score >= 0:
#            raise ValueError("score should be >= 0")
#        self._allele = allele
#        self._start = start
#        self._end = end
#        self._score = score
#
#    @property
#    def allele(self):
#        return self._allele
#
#    @property
#    def start(self):
#        return self._start
#
#    @property
#    def end(self):
#        return self._end
#
#    @property
#    def score(self):
#        return self._score
#
#    def __eq__(self, other):
#        if isinstance(other, IntervalAnnotation):
#            return self._allele == other._allele and \
#                    self._start == other._start and \
#                    self._end == other._end and \
#                    self._score == other._score
#        else:
#            raise NotImplementedError()

class ConsensusIntervalAnnotation(IntervalAnnotation):
    """Holds consensus annotation for an interval on a sequence.

    The consensus allele is decided by the AnnotationInterval with highest
    count, then by score (assuming higher score is better). The consensus
    interval is decided similarly, but is based on only the consensus allele.
    """
    def __init__(self, other = None, count = 1):
        super(ConsensusIntervalAnnotation, self).__init__(None, 0, 0, 0)
        self._alleles = {}
        self._table = {}
        if not other is None:
            self.add(other, count)

    def _add_interval_annotation(self, ia, count):
        assert count >= 0
        self._alleles[ia.allele.name] = ia.allele

        # update allele frequency & score
        afs = self._table.setdefault(ia.allele.name, [0, 0, {}])
        afs[0] += count
        afs[1] = max(afs[1], ia.score)

        # update interval frequency & score for the allele
        ifs = afs[2].setdefault((ia.start, ia.end), [0, 0])
        ifs[0] += count
        ifs[1] = max(ifs[1], ia.score)

    def _add_consensus_interval_annotation(self, cia):
        self._alleles.update(cia._alleles)
        for name, other_afs in cia._table.iteritems():

            # Update allele frequency & score
            afs = self._table.setdefault(name, [0, 0, {}])
            afs[0] += other_afs[0]
            afs[1] = max(afs[1], other_afs[1])

            # Update interval frequency & score for the allele
            for other_interval, other_ifs in other_afs[2].iteritems():
                ifs = afs[2].setdefault(other_interval, [0, 0])
                ifs[0] += other_ifs[0]
                ifs[1] = max(ifs[1], other_ifs[1])

    def add(self, other, count = 1):
        """Add IntervalAnnotation-like object and update consensus accordingly.

        :other: IntervalAnnotation-like object
        :count: frequency of 'other', should be >= 1. This paramter is ignored
        if 'other' is a ConsensusIntervalAnnotation object.
        """
        if other._allele is None:
            return
        if not self._allele is None:
            if other._allele.gene != self._allele.gene:
                raise ValueError(
                        "trying to add %s gene allele to %s gene"%(
                            other._allele.gene, self._allele.gene))
            if other._allele.region != self._allele.region:
                raise ValueError(
                        "trying to add %s region allele to %s region"%(
                            other._allele.region, self._allele.region))
        if isinstance(other, ConsensusIntervalAnnotation):
            self._add_consensus_interval_annotation(other)
        else:
            self._add_interval_annotation(other, count)
        self._update()
        return self

    def _update(self):
        # Get the allele with highest frequency and score
        name, afs = max(self._table.iteritems(), key = itemgetter(1))

        # Get interval with highest frequency and score
        (start, end), ifs = max(afs[2].iteritems(), key = itemgetter(1))

        self._allele = self._alleles[name]
        self._start = start
        self._end = end
        self._score = afs[1]
