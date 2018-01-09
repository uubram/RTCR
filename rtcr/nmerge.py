#!/usr/bin/env python
# Description: merges clones with N in the sequence to the nearest clone
# (Hamming-distance-wise).

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

from itertools import chain, groupby

from clone import SearchableCloneSet

def run_nmerge_on_bin(cloneset):
    if isinstance(cloneset, SearchableCloneSet):
        scls = cloneset
    else:
        scls = SearchableCloneSet(cloneset)
    del cloneset

    if scls.count == 0:
        return scls

    seqlen = len(next(iter(scls)).seq)

    logger.info("NMerge start (seqlen: %s)"%seqlen)
    for clone in sorted(scls, key = \
            lambda clone:(min(clone.qual), clone.count)):
        if not 'N' in clone.seq:
            continue
        if not clone in scls:
            continue

        for hd, nb in sorted(scls.neighbors(clone,
            maxhd = clone.seq.count('N')), key = \
                    lambda (hd, nb):(hd, min(nb.qual), nb.count)):
            if not nb in scls:
                continue
            if 'N' in nb.seq:
                continue
           
            mc_diff, merged_clone = scls.premerge(clone, nb)
            if not merged_clone is None:
                scls.enact_merge(clone, nb, merged_clone)
                logger.debug("merged %s to %s"%(clone, nb))
                break
    logger.info("NMerge end (seqlen: %s)"%seqlen)
    return scls

def main():
    import argparse
    import cPickle as pickle
    from fileio import LegacyCloneSetFormat
    from clone import CloneSet
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required = True,
            help = "Python object file with a CloneSet object")
    parser.add_argument("-o", required = True, help = "output tsv file.")
    args = parser.parse_args()

    input_fn = args.i
    output_fn = args.o

    logging.basicConfig(level = logging.DEBUG)
   
    print "Getting cloneset"
    with open(input_fn, 'r') as infile:
        cloneset = pickle.load(infile)

    by_seqlen = lambda clone:len(clone.seq)
    results = []
    for seqlen, clones in groupby(sorted(cloneset, key = by_seqlen),
            by_seqlen):
        results.append(run_nmerge_on_bin(CloneSet(clones)))
    
    cloneset = CloneSet(chain.from_iterable(results))

    print "Writing cloneset"
    LegacyCloneSetFormat.records_out(open(output_fn, 'w'), cloneset)
    print "\nClones left: %s"%len(cloneset)

if __name__ == "__main__":
    main()
