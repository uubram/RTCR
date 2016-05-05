#!/usr/bin/env python
# Description: Attempts to resolve sequencing substitution errors.

import logging
logger = logging.getLogger(__name__)

from math import ceil, log10, sqrt, isnan
from collections import namedtuple
from itertools import izip, ifilter
from heap import Heap
from util import binom_pmf, norm_ppf
from clone import SearchableCloneSet, clone_merge_stats

MergeTask = namedtuple("MergeTask", ["Qm", "hd", "neg_max_orig_count",
"mc_diff", "clone1", "clone2", "merged_clone"])

def add_task(scls, clone2task, tasks, clone1, clone2, max_Qm, max_Q):
    if clone1 == clone2:
        return

    mc_diff, merged_clone = scls.premerge(clone1, clone2)

    if merged_clone is None:
        return

    hd, Qm, max_mut_q = clone_merge_stats(clone1, clone2)

    if Qm > max_Qm:
        return

    if max_mut_q >= max_Q:
        return
    
    task = MergeTask(
            Qm = Qm,
            hd = hd,
            neg_max_orig_count = -max(clone1.orig_count, clone2.orig_count),
            mc_diff = mc_diff,
            clone1 = clone1,
            clone2 = clone2,
            merged_clone = merged_clone)

    # Clones that are removed by the task, have to be associated with the task
    for clone in ifilter(lambda clone: clone != merged_clone,
            (clone1, clone2)): 
        if clone in clone2task and task >= clone2task[clone]:
            # a better or equally good task already exists
            return
    for clone in ifilter(lambda clone: clone != merged_clone,
            (clone1, clone2)): 
        if clone in clone2task:
            old_task = clone2task[clone]
            tasks.remove_item(old_task)
            for otc in ifilter(lambda otc: otc != old_task.merged_clone,
                    (old_task.clone1, old_task.clone2)):
                if otc in clone2task:
                    del clone2task[otc]

        if not task in tasks:
            tasks.push(task)
            clone2task[clone] = task

def run_qmerge_on_bin(cloneset, mismatch_rate, confidence, max_Q):
    if mismatch_rate == 0 or len(cloneset) == 0:
        return cloneset

    assert(mismatch_rate < 1.0)

    seq_len = len(iter(cloneset).next().seq)
    z_conf = norm_ppf(confidence)
    n = float(cloneset.base_count)
    p = float(mismatch_rate)
    max_mutation_count = n*p + z_conf*sqrt(n*p*(1-p))

    # Calculate binomial pmf for current sequence length and mismatch_rate.
    hd2prob = [ binom_pmf(hd, seq_len, mismatch_rate)
            for hd in xrange(seq_len + 1) ]
 
    hd2nes = [cloneset.sequence_count * prob for prob in hd2prob]

    max_Qm = int(ceil(10 * log10(cloneset.base_count)))
    max_hd = max([i if r >= 1.0 else float("nan") \
            for i, r in enumerate(hd2nes)])
    if isnan(max_hd) or max_hd < 1:
        logger.info("maxhd < 1, not performing QMerge (len: %s)"%seq_len)
        return cloneset

    if isinstance(cloneset, SearchableCloneSet):
        scls = cloneset
    else:
        scls = SearchableCloneSet(cloneset)
    del cloneset

    tasks = Heap()
    basemsg = lambda:"QMerge (len: %s): "%seq_len
    fullbasemsg = lambda:"QMerge (len: %s;max_hd: %s, max_Qm: %s;max_Q: %s;\
#clones: %s;#sequences: %s;#bases: %s;#mutations: %s;max_mutations: %s;\
#tasks: %s): "%(
            seq_len, max_hd, max_Qm, max_Q, scls.count, scls.sequence_count,
            scls.base_count, scls.mutation_count, max_mutation_count,
            len(tasks))
    logger.info(fullbasemsg() + "start")
    clone2task = {}
    for hd, clone1, clone2 in scls.pairs(seq_len, max_hd):
        add_task(scls, clone2task, tasks, clone1, clone2, max_Qm, max_Q)
    logger.debug(fullbasemsg() + "finished searching for pairs")

    while len(tasks) > 0 and scls.mutation_count <= max_mutation_count:
        task = tasks.pop()
        clone1, clone2 = task.clone1, task.clone2

        # Remove current task from clone2task dict to keep it in sync with
        # the tasks heap.
        if clone1 in clone2task and task == clone2task[clone1]:
            del clone2task[clone1]
        if clone2 in clone2task and task == clone2task[clone2]:
            del clone2task[clone2]

        # Check if task can be performed
        if not clone1 in scls or not clone2 in scls:
            # Task cannot be performed, need to search for an alternate
            # "best task" for one or both clones.
            for clone in (clone1, clone2):
                if merged_clone != clone and clone in scls:
                    # New tasks may need to be created
                    for hd, neighbor in scls.neighbors(clone, max_hd):
                        add_task(scls, clone2task, tasks, clone, neighbor,
                                max_Qm, max_Q)
            continue

        # Recalculate merge because clones may have changed by other merges.
        mc_diff, merged_clone = scls.premerge(clone1, clone2)

        if merged_clone is None:
            continue

        if scls.mutation_count + mc_diff > max_mutation_count:
            continue

        hd, Qm, max_mut_q = clone_merge_stats(clone1, clone2)
        if Qm > max_Qm:
            continue

        if max_mut_q >= max_Q:
            continue

        logger.debug(basemsg() + "merging: %s(%s) %s(%s)"%(
            clone1.seq, clone1.count, clone2.seq, clone2.count))

        new_clone = not merged_clone in scls
        scls.enact_merge(clone1, clone2, merged_clone)

        if new_clone:
            logger.debug(basemsg() + "new clone: %s"%(merged_clone.seq))
            for hd, neighbor in scls.neighbors(merged_clone, max_hd):
                add_task(scls, clone2task, tasks, merged_clone, neighbor,
                        max_Qm, max_Q)

    logger.info(fullbasemsg() + "end")

    return scls

def main():
    import argparse
    from rtcr_io import LegacyCloneSetFormat
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",required = True,\
            help = "tsv input file with candidates")
    parser.add_argument("-o",required = True,help = "output tsv file.")
    parser.add_argument("-m","--mismatch_rate",required = True,\
            type = float, help = "Per base mismatch rate")
    parser.add_argument("-p",type = int, help="dummy argument")
    args = parser.parse_args()

    input_fn = args.i
    output_fn = args.o
    mismatch_rate = args.mismatch_rate

    if mismatch_rate < 0.0 or mismatch_rate > 1.0:
        raise Exception("Mismatch rate should be in range [0,1]")

    print "Getting cloneset"
    cloneset = LegacyCloneSetFormat.records_in(open(input_fn, 'r'))

    max_Q = max([max(clone.qual) for clone in cloneset])
    cloneset = run_qmerge_on_bin(cloneset, mismatch_rate, confidence = .5,
            max_Q = max_Q)

    print "Writing cloneset"
    LegacyCloneSetFormat.records_out(open(output_fn, 'w'), cloneset)

    print "\nCandidates left: %s"%len(cloneset)
    print "Bases changed: %s"%cloneset.mutation_count
    print "max bases changed: %s"%(cloneset.base_count * mismatch_rate)
        
if __name__ == "__main__":
    main()
