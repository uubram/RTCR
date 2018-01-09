#!/usr/bin/env python

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

from itertools import groupby
from math import ceil, log10, isnan

from editdistance import eval as ed 

from imerge import MergeInfo
from clone import CloneSet, Clone, AnnotatedClone
from util import binom_pmf, norm_ppf

def run_lmerge(cloneset, mismatch_rate, insertion_rate, deletion_rate,
        confidence):
    if len(cloneset) == 0 or \
            (mismatch_rate == insertion_rate == deletion_rate == 0):
        return cloneset

    assert(0 <= mismatch_rate < 1.0)
    assert(0 <= insertion_rate < 1.0)
    assert(0 <= deletion_rate < 1.0)
    assert(0 < confidence < 1.0)

    max_mutation_count = float(cloneset.base_count) * mismatch_rate
    max_insertion_count = float(cloneset.base_count) * insertion_rate
    max_deletion_count = float(cloneset.base_count) * deletion_rate

    logger.info("(LMerge start: ) mismatch_rate: %(mismatch_rate)s, \
insertion_rate: %(insertion_rate)s, deletion_rate: %(deletion_rate)s\n\
max_mutation_count: %(max_mutation_count)s, \
max_insertion_count: %(max_insertion_count)s, \
max_deletion_count: %(max_deletion_count)s"%locals())

    z_conf = norm_ppf(confidence)

    insertion_count = 0
    deletion_count = 0

    oof2if_only = True # only perform out-of-frame to in-frame merges

    # Divide the cloneset into length families (i.e. clones with equal length
    # sequences)
    by_seqlen = lambda clone:len(clone.seq)
    lenfams = {}
    for seqlen, clones in groupby(sorted(cloneset, key = by_seqlen),
            by_seqlen):
        lenfams[seqlen] = CloneSet(clones)

    clone2mergeinfo = {}
    mismatch_probs = {}

    # TODO: calculate more precise upper bound on the maximum length difference
    max_lendiff = max(lenfams.iterkeys())
    logger.info("LMerge max length difference: %s"%max_lendiff)
    for lendiff in xrange(1, max_lendiff + 1):
        logger.info("LMerge: processing length difference %s"%lendiff)
        for seqlen in sorted(lenfams, reverse = True):
            lenfam = lenfams[seqlen]

            if oof2if_only and seqlen % 3 != 0:
                # skip out-of-frame length family
                logger.debug(
                        "Skipping out-of-frame length family (seqlen: %s)"%\
                                seqlen)
                continue

            del_prob = binom_pmf(lendiff, seqlen, deletion_rate)
            # note: there is one extra place for an insertion to take place
            # than there are bases in the sequences (hence the seqlen + 1).
            ins_prob = binom_pmf(lendiff, seqlen + 1, insertion_rate)

            max_indel_prob = max(ins_prob, del_prob)

            if 1.0 / lenfam.sequence_count > max_indel_prob:
                logger.debug(
                        "Expecting < 1 indels for length family (seqlen: %s)"%\
                                seqlen)
                continue
            else:
                # TODO: it might be better to calculate what the estimated
                # true clone count should be. Since then clones that had a
                # 0 count in the raw data are then enabled to absorb indel
                # sequences (though this is probably an exceedingly rare event)

                # Since we will compare against clone.orig_count values,
                # we need to calculate the minimum orig_count
                min_clone_orig_count = (1.0 / max_indel_prob) * \
                        binom_pmf(0, seqlen, mismatch_rate)
                logger.debug(
                        "(seqlen: %s; lendiff: %s) min_clone_orig_count: %s"%(\
                        seqlen, lendiff, min_clone_orig_count))

            dels_lenfam = lenfams.get(seqlen - lendiff, CloneSet())
            ins_lenfam = lenfams.get(seqlen + lendiff, CloneSet())

            if dels_lenfam.count == 0 and ins_lenfam.count == 0:
                logger.debug("lenfams %s and %s contain no sequences"%(
                    seqlen - lendiff, seqlen + lendiff))
                continue

            for clone in sorted(lenfam, key = lambda clone:-clone.orig_count):
                if clone.orig_count < min_clone_orig_count:
                    break
                if not clone in cloneset:
                    continue
                if oof2if_only and not len(clone.seq) % 3 == 0:
                    continue

                hd2prob = mismatch_probs.setdefault(seqlen,
                        [ binom_pmf(hd, seqlen, mismatch_rate) \
                                for hd in xrange(seqlen + 1) ])
                mergeinfo = clone2mergeinfo.setdefault(clone,
                        MergeInfo(clone, hd2prob, z_conf))

                # Maximum levenshtein distance
                max_ld = lendiff + mergeinfo.maxhd

                for nb_lenfam, p in ((dels_lenfam, del_prob),
                        (ins_lenfam, ins_prob)):
                    max_neighbor_count = ceil(p * clone.count)

                    # TODO: calculate banded Levenshtein distance for speedup?
                    nbs = [ (nb, ed(nb.seq, clone.seq)) for nb in nb_lenfam]
                    nbs = [ (nb, levdist) for nb,levdist in nbs \
                            if levdist <= max_ld ]

                    for nb, levdist in sorted(nbs, key = \
                            lambda (nb, levdist):(levdist, nb.count,
                                min(nb.qual))):
                        hd = levdist - lendiff
                        if hd * nb.count + cloneset.mutation_count > \
                                max_mutation_count and hd > 0:
                            logger.debug("Skipping merge: \
hd (%s) * nb.count (%s) + cloneset.mutation_count (%s) > \
max_mutation_count(%s)"%(hd, nb.count, cloneset.mutation_count,
    max_mutation_count))
                            continue
                        if nb.seq > seqlen:
                            if lendiff*nb.count + deletion_count > \
                                    max_deletion_count:
                                logger.debug("Skipping merge: \
lendiff (%s) * nb.count (%s) + deletion_count (%s) > \
max_deletion_count (%s)"%(lendiff, nb.count, deletion_count,
    max_deletion_count))
                                continue
                        else:
                            if lendiff*nb.count + insertion_count > \
                                    max_insertion_count:
                                logger.debug("Skipping merge: \
lendiff (%s) * nb.count (%s) + insertion_count (%s) > \
max_insertion_count(%s)"%(lendiff, nb.count, insertion_count,
    max_insertion_count))
                                continue

                        if nb.count > max_neighbor_count:
                            continue
                        if not nb in cloneset:
                            continue

                        # TODO: do a proper merge here, calculating pairwise
                        # sequence alignment

                        # Since the clone.count cannot be modified directly,
                        # first create a new neighbor with the same
                        # sequence as clone, and then add this new neighbor
                        # to clone to update its count that way.
                        try:
                            nb2 = type(nb)(nb.v, nb.d, nb.j, nb.c,
                                    other = (clone.seq, clone.qual, nb.count))
                        except: # apparently nb is not an AnnotatedClone
                            nb2 = type(nb)((clone.seq, clone.qual, nb.count))

                        new_clone = clone.add(nb2)
                        if not new_clone is None:
                            cloneset.remove(nb)
                            cloneset.remove(clone)
                            cloneset.add(new_clone)
                            if len(nb.seq) > clone.seq:
                                deletion_count += lendiff * nb.count
                            else:
                                insertion_count += lendiff * nb.count
                        logger.debug("merged: %s to %s"%(nb, clone))

    return cloneset

def main():
    import argparse
    import cPickle as pickle
    from fileio import LegacyCloneSetFormat
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required = True,
            help = "Python object file with a CloneSet object")
    parser.add_argument("-o", required = True, help = "output tsv file.")
    parser.add_argument("-mr", "--mismatch_rate", required = True,
            type = float, help = "Per base mismatch rate")
    parser.add_argument("-ir", "--insertion_rate", required = True,
            type = float, help = "Per base insertion rate")
    parser.add_argument("-dr", "--deletion_rate", required = True,
            type = float, help = "Per base deletion rate")
    parser.add_argument("-c", "--confidence", required = True, type=float,
            help = "probability actual size of a candidate is \
less than the calculated estimate.")
    args = parser.parse_args()

    input_fn = args.i
    output_fn = args.o
    mismatch_rate = args.mismatch_rate
    insertion_rate = args.insertion_rate
    deletion_rate = args.deletion_rate
    confidence = args.confidence

    assert(0.0 < mismatch_rate < 1.0)
    assert(0.0 < insertion_rate < 1.0)
    assert(0.0 < deletion_rate < 1.0)
    assert(0.0 < confidence < 1.0)

    logging.basicConfig(level = logging.DEBUG)
   
    print "Getting cloneset"
    with open(input_fn, 'r') as infile:
        cloneset = pickle.load(infile)
    
    cloneset = run_lmerge(cloneset, mismatch_rate, insertion_rate,
            deletion_rate, confidence)

    print "Writing cloneset"
    LegacyCloneSetFormat.records_out(open(output_fn, 'w'), cloneset)
    print "\nClones left: %s"%len(cloneset)
 
if __name__ == "__main__":
    main()
