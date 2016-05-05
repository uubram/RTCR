#!/usr/bin/env python
# Description: Attempts to resolve sequencing and PCR substitution errors by
# merging candidates based on size difference and Hamming Distance.

import logging
logger = logging.getLogger(__name__)

from math import isnan, sqrt
from collections import defaultdict
from itertools import chain
from util import binom_pmf, norm_ppf
from clone import CloneSet, SearchableCloneSet

class MergeInfo():
    """Contains information required by the IMerge algorithm to perform merges
    on a particular clone.

    This is used by IMerge algorithm to determine how many and which sequences
    in the sequence space neighborhood of a particular clone can be merged to
    the clone.
    """
    def __init__(self, clone, hd2prob, z_conf):
        p = hd2prob[0] # probability for a correct sequence
        est_clonesize = \
                (clone.orig_count + z_conf*sqrt(clone.orig_count*(1-p))) / p
        
        # Calculate number of erroneous sequences expected in the neighborhood
        # of the current clone (ignoring contributions from other clones).
        total_remaining = max(est_clonesize - clone.count, 0)

        err_frac = 1-p

        # Calculate number of erroneous sequences expected to be at every HD,
        # given the difference between original and current clonesize.
        hd2nes = [0] + [total_remaining * (prob/err_frac) \
                for prob in hd2prob[1:]]

        # Calculate first Hamming distance at which less than 1 erroneous
        # sequence is expected.
        self.maxhd = max([i if est_clonesize * prob >= 1.0 else float("nan")\
                for i, prob in enumerate(hd2prob)]) + 1

        if isnan(self.maxhd):
            self.maxhd = 0
        
        self._neighborhood = None

        self.clone = clone
        self.hd2nes = hd2nes
        # sweep: Cumulative number of unmerged erroneous sequences, 
        # counted from Hamming distance 1
        self.sweep = 0

    def neighborhood(self, scls, extend_hd):
        if self._neighborhood is None or extend_hd > self.maxhd:
            if extend_hd > self.maxhd:
                self.maxhd = extend_hd 
            self._neighborhood = set(scls.neighbors(self.clone, self.maxhd))
        return self._neighborhood

def run_imerge(cloneset, mismatch_rate, confidence):
    bins = defaultdict(lambda:CloneSet())
    for clone in cloneset:
        bins[len(clone.seq)].add(clone)
    return CloneSet(chain.from_iterable([
            run_imerge_on_bin(cloneset, mismatch_rate, confidence) \
                    for cloneset in bins.itervalues()]))

#TODO: add listener
def run_imerge_on_bin(cloneset, mismatch_rate, confidence):
    """Run IMerge algorithm.

    :param cloneset: a set of clones of equal length sequences.
    :param mismatch_rate: probability any given base in the cloneset is
    erroneous. float value with range [0,1].
    :confidence: confidence level used for clone size estimation. 
    """

    if mismatch_rate == 0 or len(cloneset) == 0:
        return cloneset

    assert(mismatch_rate < 1.0)
    assert(0.0 < confidence and confidence < 1.0)

    seq_len = len(iter(cloneset).next().seq)
    z_conf = norm_ppf(confidence)
    n = float(cloneset.base_count)
    p = float(mismatch_rate)
    max_mutation_count = n*p + z_conf*sqrt(n*p*(1-p))

    if isinstance(cloneset, SearchableCloneSet):
        scls = cloneset
    else:
        scls = SearchableCloneSet(cloneset)
    del cloneset

    # Calculate binomial pmf for current sequence length and mismatch_rate.
    hd2prob = [ binom_pmf(hd, seq_len, mismatch_rate)
            for hd in xrange(seq_len + 1) ]
    clone2mergeinfo = {clone : MergeInfo(clone, hd2prob, z_conf) \
            for clone in scls}

    cycle_maxhd = max([mergeinfo.maxhd for mergeinfo in \
            clone2mergeinfo.itervalues()])

    basemsg = lambda:"IMerge (len: %s): "%seq_len
    fullbasemsg = lambda:"IMerge (len: %s;cycle_maxhd: %s;\
#clones: %s;#sequences: %s;#bases: %s;#mutations: %s;max_mutations: %s): "%(
            seq_len, cycle_maxhd, scls.count, scls.sequence_count,
            scls.base_count, scls.mutation_count, max_mutation_count)

    if isnan(cycle_maxhd):
        logger.info(fullbasemsg + "stopping because cycle_maxhd is 0"%\
                seq_len)
        return scls
   
    logger.info(fullbasemsg() + "start")
    for cycle_hd in xrange(1, cycle_maxhd + 1):
        logger.debug(basemsg() + "cycle %s"%cycle_hd)
        for clone in sorted(scls, key = lambda clone : \
                (-clone.count, clone.seq)):
            if not clone in scls:
                logger.debug(basemsg() + "%s not in scls"%clone.seq)
                continue

            mergeinfo = clone2mergeinfo[clone]

            if cycle_hd > mergeinfo.maxhd:
                logger.debug(basemsg() + "%s (maxhd: %s; cycle_hd: %s)"%(
                    clone.seq, mergeinfo.maxhd, cycle_hd))
                continue

            mergeinfo.sweep += mergeinfo.hd2nes[cycle_hd]

            logger.debug(basemsg() + "seq: %s (%s; sweep: %s)"%(clone.seq,
                clone.count, mergeinfo.sweep))

            for hd, neighbor in sorted(mergeinfo.neighborhood(scls, cycle_hd),
                    key = lambda (hd, neighbor): (hd, neighbor.count)):
                if clone == neighbor or neighbor not in scls:
                    continue
                if hd > cycle_hd:
                    break
                if mergeinfo.sweep >= neighbor.count:
                    mc_diff, merged_clone = scls.premerge(clone, neighbor)

                    if merged_clone is None:
                        continue

                    if scls.mutation_count + mc_diff > max_mutation_count:
                        logger.debug(fullbasemsg() + "mc_diff = %s,\
clone.count = %s, neighbor.count = %s"%(mc_diff, clone.count, neighbor.count))
                        continue

                    logger.debug(basemsg() + "merging: %s(%s) %s(%s)"%(
                        clone.seq, clone.count, neighbor.seq, neighbor.count))
                    scls.enact_merge(clone, neighbor, merged_clone)

                    mergeinfo.sweep -= neighbor.count
                    del clone2mergeinfo[neighbor]
                    del clone2mergeinfo[clone]
                    mergeinfo.clone = merged_clone
                    clone2mergeinfo[merged_clone] = mergeinfo
                    clone = merged_clone
    logger.info(fullbasemsg() + "end")
    return scls

def main():
    import argparse
    from fileio import LegacyCloneSetFormat
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",  required = True,
            help = "tsv input file with candidates")
    parser.add_argument("-o",required = True,help = "output tsv file.")
    parser.add_argument("-mr","--mismatch_rate", required = True,
            type = float, help = "Per base mismatch rate")
    parser.add_argument("-c","--confidence",required = True,type=float,
            help = "probability actual size of a candidate is \
less than the calculated estimate.")
    args = parser.parse_args()

    input_fn = args.i
    output_fn = args.o
    mismatch_rate = args.mismatch_rate
    confidence = args.confidence

    assert(0.0 < mismatch_rate < 1.0)
    assert(0.0 < confidence < 1.0)
   
    print "Getting cloneset"
    cloneset = LegacyCloneSetFormat.records_in(open(input_fn, 'r'))

    cloneset = run_imerge(cloneset, mismatch_rate, confidence)

    print "Writing cloneset"
    LegacyCloneSetFormat.records_out(open(output_fn, 'w'), cloneset)

    print "\nCandidates left: %s"%len(cloneset)
    print "Bases changed: %s"%cloneset.mutation_count
    print "max bases changed: %s"%(cloneset.base_count * mismatch_rate)
 
if __name__ == "__main__":
    main()
