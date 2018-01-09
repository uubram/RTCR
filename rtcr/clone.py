#!/usr/bin/env python
#
# Description: provides datastructures to hold (B/T)cR clones.

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

from collections import deque, Iterable, namedtuple
from itertools import izip
from allele import AlleleContainer
from seq import merge_stats, cigar_intervals, cigar_rpos2qpos, \
        ConsensusQSequence, IntervalAnnotation, ConsensusIntervalAnnotation
from trie import Trie

def clone_merge_stats(clone1, clone2):
    return merge_stats(clone1._cqs, clone2._cqs, clone1.refpos - clone2.refpos)

class Clone(object):
    """Immutable set of identical sequences.

    A clone is uniquely identified by its (consensus) sequence. For example,
    two Clone objects with different refpos, orig_count, and count are still
    considered the same clone if their sequences are the same. Derived classes
    are *not* expected to adhere to this and can define uniqueness of a clone
    on more than sequence alone.
    """

    __slots__ = ['_cqs', '_seq', '_qual', '_count', '_base_count',
            '_mutation_count', '_refpos', '_orig_count']

    def __getstate__(self):
        return self._cqs, self._seq, self._qual, self._count, \
                self._base_count, self._mutation_count, self._refpos, \
                self._orig_count
    def __setstate__(self, state):
        self._cqs, self._seq, self._qual, self._count, \
                self._base_count, self._mutation_count, self._refpos, \
                self._orig_count = state

    def __init__(self, other, refpos = 0, orig_count = None):
        """
        :other: a Clone, ConsensusQSequence, QSequence-like, or list/tuple
        object containing in order seq, qual, and count values.
        :refpos: reference position, this is used to align clones when adding
        them together.
        :orig_count: original count of the clone when it was first created.
        """
        if isinstance(other, Clone):
            self._cqs = ConsensusQSequence(other._cqs)
            if orig_count is None:
                orig_count = other._orig_count
        else:
            self._cqs = ConsensusQSequence(other)
            if orig_count is None:
                orig_count = self._cqs.count

        assert orig_count >= 0
        assert orig_count <= self._cqs.count

        self._sync_with_cqs()
        self._refpos = refpos
        self._orig_count = orig_count

    def _sync_with_cqs(self):
        self._seq = self._cqs.seq
        self._qual = self._cqs.qual
        self._count = self._cqs.count
        self._base_count = self._cqs.base_count
        self._mutation_count = self._cqs.mutation_count

    @classmethod
    def from_pfm(cls):
        raise NotImplementedError()

    @property
    def seq(self):
        return self._seq

    @property
    def qual(self):
        return self._qual

    @property
    def count(self):
        return self._count

    @property
    def base_count(self):
        return self._base_count

    @property
    def mutation_count(self):
        return self._mutation_count

    @property
    def orig_count(self):
        """Return the number of sequences."""
        return self._orig_count

    @property
    def refpos(self):
        return self._refpos

    def add(self, other):
        """Return new clone resulting from adding other to self."""

        clone = self.copy()
        clone._cqs.add(other._cqs, offset = self._refpos - other._refpos)
        clone._refpos = max(self._refpos, other._refpos)
        clone._orig_count = max(self.orig_count, other.orig_count)
        clone._sync_with_cqs()
        return clone
    
    def copy(self):
        """Create copy deep enough that modifying attributes of the copy does
        not change the original.
        """

        cls = self.__class__
        cpy = cls.__new__(cls)
        for name in cpy.__slots__:
            if name != "_cqs":
                setattr(cpy, name, getattr(self, name))
        cpy._cqs = ConsensusQSequence(self._cqs)
        return cpy

    def __len__(self):
        return len(self._seq)

    def __str__(self):
        return self._seq

    def __eq__(self, other):
        if isinstance(other, Clone):
            return self._seq == other._seq
        else:
            raise NotImplementedError()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self._seq)

class AnnotatedClone(Clone):
    """Represents immutable nucleotide sequence of a single (B/T)-cell receptor
    chain.

    This class uniquely identifies a clone the same way its parent class does.
    There are no restrictions placed on the annotation other than that a V and
    J annotation be provided, D and C are allowed to lack annotation.
    """
    __slots__ = ['_v', '_d', '_j', '_c']

    def __getstate__(self):
        super_state = super(AnnotatedClone, self).__getstate__()
        return (self._v, self._d, self._j, self._c, super_state)

    def __setstate__(self, state):
        self._v, self._d, self._j, self._c, super_state = state
        super(AnnotatedClone, self).__setstate__(super_state)

    def __init__(self, v, d, j, c, *args, **kwargs):
        """
        :v: IntervalAnnotation-like object for V-region of the clone sequence.
        :d: similar, but for D-region
        :j: J-region
        :c: C-region
        """
        super(AnnotatedClone, self).__init__(*args, **kwargs)
        if v is None or v.allele is None or v.allele.region != "V":
            raise ValueError("v param is not a V-region")
        if not d is None and not d.allele is None and d.allele.region != "D":
            raise ValueError("d param is not a D-region")
        if j is None or j.allele is None or j.allele.region != "J":
            raise ValueError("j param is not a J-region")
        if not c is None and not c.allele is None and  c.allele.region != "C":
            raise ValueError("c param is not a C-region")
        self._v = v
        self._d = IntervalAnnotation(None, 0, 0, 0) if d is None else d
        self._j = j
        self._c = IntervalAnnotation(None, 0, 0, 0) if c is None else c

    def add(self, other):
        if not isinstance(other, self.__class__):
            raise ValueError("not an instance of %s"%self.__class__)

        clone = super(AnnotatedClone, self).add(other)

        if other.count > self.count or \
                (other.count == self.count and \
                other.v.score + other.d.score + other.j.score + \
                other.c.score > self.v.score + self.d.score +self.j.score + \
                self.c.score):
            annot = other
        else:
            annot = self

        clone._v = annot._v
        clone._d = annot._d
        clone._j = annot._j
        clone._c = annot._c

        return clone
    
    @property
    def v(self):
        return self._v

    @property
    def d(self):
        return self._d

    @property
    def j(self):
        return self._j

    @property
    def c(self):
        return self._c

class AnnotatedCloneDistinctAllele(AnnotatedClone):
    """Represents immutable nucleotide sequence of a single (B/T)-cell receptor
    chain.

    This class uniquely identifies a clone by sequence and V(D)J annotation.
    """

    @staticmethod
    def matches(x, y):
        if x is None or y is None:
                return True
        else:
            return x.allele == y.allele

    def add(self, other):
        if self.matches(other.v, self.v) and \
                self.matches(other.d, self.d) and \
                self.matches(other.j, self.j):
            return super(AnnotatedCloneDistinctAllele, self).add(other)
        else:
            return None

    def __str__(self):
        vid = self.v.allele.name
        did = None if self.d.allele is None else self.d.allele.name
        jid = self.j.allele.name
        cid = None if self.c.allele is None else self.c.allele.name
        return "%s(%s,%s,%s,%s,%s)"%(self.__class__.__name__, vid, did, jid,
                cid, self.seq)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self._seq == other._seq and \
                    self.matches(self.v, other.v) and \
                    self.matches(self.d, other.d) and \
                    self.matches(self.j, other.j) and \
                    self.matches(self.c, other.c)
        else:
            raise NotImplementedError()

class ConsensusAnnotatedClone(Clone):
    """Represents immutable nucleotide sequence of a single (B/T)-cell receptor
    chain.

    This class uniquely identifies a clone the same way its parent class does.
    There are no restrictions placed on the annotation other than that a V and
    J annotation be provided, D and C are allowed to lack annotation.
    """
    def __init__(self, v, d, j, c, *args, **kwargs):
        """
        :v: IntervalAnnotation-like object for V-region of the clone sequence.
        :d: similar, but for D-region
        :j: J-region
        :c: C-region
        """
        super(AnnotatedClone, self).__init__(*args, **kwargs)
        if v is None or v.allele is None or v.allele.region != "V":
            raise ValueError("v param is not a V-region")
        if not d is None and not d.allele is None and d.allele.region != "D":
            raise ValueError("d param is not a D-region")
        if j is None or j.allele is None or j.allele.region != "J":
            raise ValueError("j param is not a J-region")
        if not c is None and not c.allele is None and  c.allele.region != "C":
            raise ValueError("c param is not a C-region")
        self._v = ConsensusIntervalAnnotation(v, self.count)
        self._d = ConsensusIntervalAnnotation(d, self.count)
        self._j = ConsensusIntervalAnnotation(j, self.count)
        self._c = ConsensusIntervalAnnotation(c, self.count)

    def add(self, other):
        if not isinstance(other, AnnotatedClone):
            raise ValueError("not an instance of AnnotatedClone")

        return AnnotatedClone(
                ConsensusIntervalAnnotation(self._v).add(other._v),
                ConsensusIntervalAnnotation(self._d).add(other._d),
                ConsensusIntervalAnnotation(self._j).add(other._j),
                ConsensusIntervalAnnotation(self._c).add(other._c),
                other = super(AnnotatedClone, self).add(other))

    @property
    def v(self):
        return self._v

    @property
    def d(self):
        return self._d

    @property
    def j(self):
        return self._j

    @property
    def c(self):
        return self._c

class CloneSet(object):
    """A set of Clone objects."""
    def __init__(self, iterable = None):
        self._clones = dict()
        self._base_count = 0
        self._mutation_count = 0
        self._sequence_count = 0

        if not iterable is None:
            for it in iterable:
                self.add(it)

    @property
    def count(self):
        """Return the number of clones."""
        return len(self._clones)

    @property
    def base_count(self):
        """Return the number of bases."""
        return self._base_count

    @property
    def mutation_count(self):
        """Return the number of mutated bases."""
        return self._mutation_count

    @property
    def sequence_count(self):
        """Return the number of sequences."""
        return self._sequence_count

    def _add_clone(self, other, merge):
        # NOTE: do not use isinstance here as it can interfere with unpickling.
        # Pickling creates identically named classes (if the pickling is not
        # done in a fresh python instance like in multiprocessing), so that the
        # unpickled objects according to isinstance are still different
        # non-unpickled objects with identical classname.
        if not hasattr(other, "orig_count"):
            raise ValueError("other (%r) is not a Clone(-like) object"%other)
        if other in self._clones:
            if merge:
                clone = self._clones[other]
                self.remove(clone)
                orig_count = clone.orig_count + other.orig_count
                other = clone.add(other)
                # TODO: do not update orig_count by accessing a private member!
                other._orig_count = orig_count
                if other is None:
                    raise Exception(
                    "Failed to merge supposedly compatible clones")
                self.add(other, merge = False)
                return
            else:
                raise KeyError("Clone (%r) already exists."%other)

        self._clones[other] = other
        self._base_count += other.base_count
        self._mutation_count += other.mutation_count
        self._sequence_count += other.count

    def add(self, other, merge = False):
        """Add clone or set of clones.
        
        :merge: true, if a clone is considered to already exist, merge it.
        false, if a clone is considered to already exist, throw a KeyError. In
        principle, this should only be set to true at the stage of building
        clones from reads. *Importantly*, the orig_count field will be the sum
        of the orig_count field of any pre-existing clone and the one to be
        added.
        """
        if isinstance(other, Iterable):
            for clone in other:
                self._add_clone(clone, merge)
        else:
            self._add_clone(other, merge)

    def remove(self, clone):
        # the argument clone may not actually be the one in the cloneset,
        # despite their hashes being identical.
        clone = self._clones[clone]
        del self._clones[clone]
        self._base_count -= clone.base_count
        self._mutation_count -= clone.mutation_count
        self._sequence_count -= clone.count

    def premerge(self, a, b):
        # The argument clones may not actually be the ones in the cloneset,
        # despite their hashes being identical.
        a = self._clones[a]
        b = self._clones[b]

        assert a != b

        c = a.add(b)

        if c is None:
            return 0, None

        mc = a.mutation_count + b.mutation_count

        if c != a and c != b and c in self:
            mc += self._clones[c].mutation_count
            c = self._clones[c].add(c)

        return c.mutation_count - mc, c

    def enact_merge(self, a, b, c):
        self.remove(a)
        self.remove(b)

        if c in self:
            self.remove(c)
        self.add(c)

    def merge(self, a, b):
        """Merge clone b to clone a and replace the clones with the merged
        clone.
        
        A tuple will be returned with the change in mutation count and the
        merged clone. If the merge fails, the tuple will contain 0 and None
        values for mutation count and merged clone respectively.
        """
        diff, c = self.premerge(a,b)
        if not c is None:
            self.enact_merge(a, b, c)

        return diff, c

    def __len__(self):
        return self.count

    def __iter__(self):
        return iter(self._clones)

    def __contains__(self, item):
        return item in self._clones

def build_clone(ref, v_rec, j_rec, include_cysphe,
        classname = "AnnotatedClone"):
    """Build clone from SAMRecord-like objects. Returns None if clone cannot be
    built due to unmapped alignment or because V and J alignments are in the
    opposite direction.

    :ref: AlleleContainer-like object containing the germline reference
    sequences that are referenced by the parameters with SAMRecord-like
    objects.
    :v_rec: SAMRecord-like object containing alignment of V segment
    :j_rec: SAMRecord-like object containing alignment of J segment
    :include_cysphe: include codons of the conserved Cys and Phe residues when
    extracting the CDR3 region from the alignments.
    """
    assert v_rec.QNAME == j_rec.QNAME
    if v_rec.FLAG & 4 or j_rec.FLAG & 4: # segment unmapped
        return None

    # assuming J segment was aligned on the read orientation resulting from
    # V segment alignment, it means that if the j_rec has the reverse
    # complement bit of the SAM flag set, that V and J alignments are not
    # in the same orientation
    if j_rec.FLAG & 16: # reverse complement
        return None

    vid = v_rec.RNAME
    jid = j_rec.RNAME

    assert ref[vid].region == "V"
    assert ref[jid].region == "J"

    # Get alignment intervals on the sequence
    j_offset = len(v_rec.SEQ) - len(j_rec.SEQ) 
    vas, vae, ref_vas, ref_vae = cigar_intervals(v_rec.CIGAR,
            v_rec.POS - 1)
    jas, jae, ref_jas, ref_jae = cigar_intervals(j_rec.CIGAR,
            j_rec.POS - 1)
    jas += j_offset
    jae += j_offset

    # Get junction from the sequence
    if include_cysphe:
        v_refpos_offset = -3
        j_refpos_offset = 3
    else:
        v_refpos_offset = 0
        j_refpos_offset = 0
    jncs = cigar_rpos2qpos(v_rec.CIGAR, v_rec.POS - 1,
            ref[vid].refpos + v_refpos_offset)
    j_rec_jnce = cigar_rpos2qpos(j_rec.CIGAR, j_rec.POS - 1,
            ref[jid].refpos + j_refpos_offset)

    if jncs is None or j_rec_jnce is None:
        return None

    jnce = j_offset + j_rec_jnce

    seq = v_rec.SEQ
    if jncs < 0 or jncs >= len(seq) or jnce < 0 or jnce > len(seq) or \
            jnce - jncs < 0:
        return None

    nt_jnc = seq[jncs:jnce]
    q_jnc = [ord(qch) - 33 for qch in v_rec.QUAL[jncs:jnce]]

    # Do not accept empty junction sequence
    if len(nt_jnc) < 1:
        return None

    v_annot = IntervalAnnotation(
        allele = ref[vid],
        start = 0,
        end = vae - jncs,
        score = v_rec.MAPQ)

    j_annot = IntervalAnnotation(
        allele = ref[jid],
        start = jas - jncs,
        end = jnce - jncs,
        score = j_rec.MAPQ)

    cloneclass = globals()[classname]
    assert issubclass(cloneclass, Clone)

    if issubclass(cloneclass, AnnotatedClone):
        return cloneclass(v = v_annot, d = None, j = j_annot, c = None,
                other = (nt_jnc, q_jnc, 1))
    else:
        return Clone((nt_jnc, q_jnc, 1))

class SearchableCloneSet(CloneSet):
    """Enables subsetting a cloneset on similarity to a given sequence or
    clone.
    """

    def __init__(self, iterable):
        self._trie = Trie()
        super(SearchableCloneSet, self).__init__(iterable)

    def _add_to_trie(self, clone):
        value = self._trie.setdefault(clone.seq, [])
        # Depending on type of clone, multiple clones might have the same
        # sequence, so node values should be a list of clones.
        value.append(clone)

    def add(self, clone):
        if clone in self:
            return
        super(SearchableCloneSet, self).add(clone)
        self._add_to_trie(clone)

    def remove(self, clone):
        if not clone in self:
            raise Exception("Cannot remove non-existent clone.")
        super(SearchableCloneSet, self).remove(clone)

        value = self._trie[clone.seq]
        value.remove(clone)
        if len(value) == 0:
            del self._trie[clone.seq]

    def nearest_neighbors(self, x, maxhd = None):
        """Get nearest neighboring clones Hamming-distance-wise.

        Yields tuples with Hamming distance and clone.

        :x: a sequence or clone
        :maxhd: maximum Hamming distance at which clones are considered
        neighbors. 
        """
        if hasattr(x, "seq"):
            seq = x.seq
        else:
            seq = x

        for hd, seq, value in self._trie.get_nearest_variants(seq, maxhd):
            if value:
                for clone in value:
                    yield hd, clone

    def neighbors(self, x, maxhd):
        """Get all neighboring clones within maxhd of a sequence or clone.

        Yields tuples with Hamming distance and clone.

        :x: a sequence or clone
        :maxhd: maximum Hamming distance at which clones are considered
        neighbors.
        """
        if hasattr(x, "seq"):
            seq = x.seq
        else:
            seq = x

        for hd, seq, value in self._trie.neighbors(seq, maxhd):
            if value:
                for clone in value:
                    yield hd, clone

    def pairs(self, seqlen, maxhd):
        """Generator function for iterating over pairs of clones within maxhd
        of each other.
        
        Yields tuples with Hamming distance and two clones.

        :seqlen: sequence length of the clones
        :maxhd: maximum Hamming distance to consider
        """
        for hd, s1, value1, s2, value2 in self._trie.pairs(keylen = seqlen,
                maxhd = maxhd):
            for clone1 in value1:
                for clone2 in value2:
                    yield hd, clone1, clone2
