#!/usr/bin/env python
# Description: merges clones with N in the sequence to the nearest clone
# (Hamming-distance-wise).

import logging
logger = logging.getLogger(__name__)

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
    logger.info("NMerge end (seqlen: %s)"%seqlen)
    return scls

if __name__ == "__main__":
    main()
