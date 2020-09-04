import logging
logger = logging.getLogger(__name__)

import os.path
import threading
import gc
import math
import gzip
from time import time, sleep
from collections import namedtuple
from itertools import izip, groupby, chain
import cPickle as pickle
import sys

from align import get_vj_alignments
from fileio import zopen, SAMFormat, LegacyCloneSetFormat, filesize
from clone import CloneSet, build_clone
from seq import get_error_stats
from util import ConnectedConsumerPool, clone2AIRRDict, nt2aa
from qmerge import run_qmerge_on_bin
from imerge import run_imerge_on_bin
from lmerge import run_lmerge
from nmerge import run_nmerge_on_bin

def p2q(p):
    return -10*math.log10(p)

def q2p(q):
    return math.pow(10, float(-q)/10)

def plot_Q_mm_stats(Q_mm_stats, fn):
    """Plots base quality scores vs observed(/effective) quality scores and
    saves it to a file.
    """
    try:
        import matplotlib
        matplotlib.use("PDF")
        import matplotlib.pyplot as plt
    except Exception as e:
        logger.warning("Unable to import matplotlib. Exception: %s"%e)
        return

    x = {"V":[], "J":[]} 
    y = {"V":[], "J":[]}
    vcolor = "green"
    jcolor = "blue"
    for region, color in (("V", vcolor), ("J", jcolor)):
        Q_mm = Q_mm_stats[region]["Q_mm"]
        Q_n = Q_mm_stats[region]["Q_n"]
        for q, (mm, n) in enumerate(zip(Q_mm, Q_n)):
            if mm > 0 and n > 0:
                mmr = float(mm) / n
                q_obs = p2q(mmr)
                x[region].append(q)
                y[region].append(q_obs)
    v_scatter = plt.scatter(x["V"], y["V"], c = vcolor)
    j_scatter = plt.scatter(x["J"], y["J"], c = jcolor)
    plt.plot([0,41], [0, 41], color = "black")
    plt.xlim((0,41))
    plt.ylim((0,41))
    plt.xlabel("Base quality score (Phred)")
    plt.ylabel("Observed quality (Phred)")
    plt.legend((v_scatter, j_scatter), ("V-region", "J-region"),
            scatterpoints = 1,
            loc = "upper left")
    plt.savefig(fn, bbox_inches = "tight")

def run_ec_on_bin(*args):
    cloneset, mismatch_rate, confidence, max_Q = args
    try:
        gc.disable()
        seqlen = len(next(iter(cloneset)).seq)
        logger.info("Starting QMerge on cloneset: seqlen: %s, \
#clones: %s, #sequences: %s, #bases: %s, #mutation_count: %s"%(seqlen,
len(cloneset), cloneset.sequence_count, cloneset.base_count,
cloneset.mutation_count))

        cloneset = run_qmerge_on_bin(cloneset, mismatch_rate, confidence,
                max_Q)
        logger.info("Starting IMerge on cloneset: seqlen: %s, \
#clones: %s, #sequences: %s, #bases: %s, #mutation_count: %s"%(seqlen,
len(cloneset), cloneset.sequence_count, cloneset.base_count,
cloneset.mutation_count))
        cloneset = run_imerge_on_bin(cloneset, mismatch_rate, confidence)
    finally:
        gc.enable()
    return cloneset, mismatch_rate

def wrapper_run_nmerge_on_bin(*args):
    return run_nmerge_on_bin(args[0])

class Pipeline(threading.Thread):
    def _save_cloneset(self, cloneset, name):
        self._listener.notify("Saving cloneset (\"%s\") to file"%name)

        #LegacyCloneSetFormat.records_out(zopen(name + ".tsv", 'w'),
        #        cloneset)

        with open(name + ".dat",'w') as f:
            pickle.dump(cloneset, f, protocol = -1)

    def run_pool(self, pool, desc):
        total_task_count = pool.task_count
        try:
            pool.start()

            self._listener.notify(("PROGRESSBAR", desc, "start"))
            while pool.task_count > 0 and not self.stopped():
                if pool.error():
                    raise Exception('Error occurred in worker pool')
                sleep(self._update_interval)
                frac = float(total_task_count - pool.task_count) / \
                        total_task_count
                self._listener.notify(("PROGRESSBAR", desc, frac))

            if self.stopped():
                raise Exception('Pipeline stopped')

            pool.join()

            self._listener.notify(("PROGRESSBAR", desc, "end"))
        except:
            pool.terminate()
            raise
        
    def __init__(self, ref, reads, phred_encoding,
            cmd_build_index,
            args_build_index,
            cmd_align,
            args_align_v,
            args_align_j,
            alignments_fn,
            alignment_stats_fn,
            Q_mm_stats_fn,
            Q_mm_stats_plot_fn,
            output_fn,
            output_not_ok_fn,
            clone_classname,
            confidence,
            min_seqlen,
            min_phred_threshold,
            n_threads,
            update_interval,
            listener):
        """
        Run Recover T Cell Receptor (RTCR) pipeline.

        Performs the following steps, starting from raw reads:
        1) aligns V and J reference sequences to every read
        2) generates clones and simultaneously estimates error rates from the
        alignments.
        3) QMerge
        4) IMerge
        5) LMerge
        6) NMerge

        .. NOTE::
           filenames that are set to None or empty string will not be written.
        
        :ref: AlleleContainer-like, should contain germline reference sequences
        :reads: handle to FastQ file
        :phred_encoding: encoding used for Phred quality scores of the reads.
        :cmd_build_index: command that will be run to build an index of the
        reference sequences. If empty "" or None, this command will not be run.
        :args_build_index: string, arguments that will be provided to
        cmd_build_index
        :cmd_align: string, command that will be run to start the aligner
        :args_align_v: string, arguments that will be provided to cmd_align to
        align V sequences to the reads.
        :args_align_j: string, arguments that will be provided to cmd_align to
        align J sequences to the reads.
        :alignments_fn: name of file to output alignments to.
        :alignment_stats_fn: name of csv file to output alignment statistics
        to, such as number of indels and mismatches.
        :Q_mm_stats_fn: name of csv file to output per Phred score mismatch
        statistics.
        :Q_mm_stats_plot_fn: name of file to output plot of observed quality vs
        base quality scores.
        :output_fn: name of file to output final clones to.
        :output_not_ok_fn: name of file to output final discarded clones to.
        :clone_classname: string, classname of type of clone to build
        :confidence: right confidence limit (between 0 and 1), this is used to
        calculate how many bases are allowed to be mutated (both QMerge and
        IMerge), and how large clones may grow (IMerge).
        :min_seqlen: minimum sequence length for a clone
        :min_phred_threshold: minimum base quality in a clone for it to be
        reported in the final results.
        :n_threads: maximum number of threads/processes the pipeline is allowed
        to use.
        :update_interval: minimum time (in seconds) between notifications to
        the listener.
        """
        super(Pipeline, self).__init__()
        self._stop = threading.Event()
        self._listener = listener
        self._ref = ref
        self._reads = reads
        self._phred_encoding = phred_encoding
        self._cmd_build_index = cmd_build_index
        self._args_build_index = args_build_index
        self._cmd_align = cmd_align
        self._args_align_v = args_align_v
        self._args_align_j = args_align_j
        self._alignments_fn = alignments_fn
        self._alignment_stats_fn = alignment_stats_fn
        self._Q_mm_stats_fn = Q_mm_stats_fn
        self._Q_mm_stats_plot_fn = Q_mm_stats_plot_fn
        self._output_fn = output_fn
        self._output_not_ok_fn = output_not_ok_fn
        self._clone_classname = clone_classname
        self._confidence = confidence
        self._min_seqlen = min_seqlen
        self._min_phred_threshold = min_phred_threshold
        self._n_threads = n_threads
        self._update_interval = update_interval

    def run(self):
        if self._alignments_fn is None or self._alignments_fn == "":
            output_alignments = False
        else:
            output_alignments = True

        if output_alignments and os.path.isfile(self._alignments_fn):
            logger.info("SKIPPING creation of %s"%self._alignments_fn)
            output_alignments = False
            alignment_file = zopen(self._alignments_fn, 'r')
            vj_recs = SAMFormat.records_in(alignment_file)
            # Get two (rows/)alignments at a time from vj_recs
            alns = ((rec, next(vj_recs)) for rec in vj_recs)
            self._listener.notify("Reading alignments from %s"%
                    self._alignments_fn)
        else:
            alns = get_vj_alignments(self._ref, self._reads,
                    self._cmd_build_index,
                    self._args_build_index,
                    self._cmd_align,
                    self._args_align_v,
                    self._args_align_j,
                    phred_encoding = self._phred_encoding,
                    n_threads = self._n_threads)
            self._listener.notify("Aligning reference sequences to reads")

        # Keep track of the quality scores of the bases that went into the
        # sequences of the clones.
        Q_counts = {}

        # Build clones and use alignments to count mismatches and indels
        cs = CloneSet()
        alnstats = {"V":{}, "J":{}}
        v_refpos_offset = -3
        j_refpos_offset = 3

        try:
            if output_alignments:
                out = zopen(self._alignments_fn, 'w')

            if output_alignments:
                infile = self._reads
            else:
                infile = alignment_file

            prev_time = time()
            if isinstance(infile, gzip.GzipFile):
                infile_size = os.path.getsize(infile.name)
                infile_pos = infile.fileobj.tell
            else:
                infile_size = filesize(infile)
                infile_pos = infile.tell

            self._listener.notify(("PROGRESSBAR", "Alignments", "start"))

            for v_rec, j_rec in alns:
                if self.stopped():
                    logger.warning("Pipeline stopped")
                    return

                if time() - prev_time >= self._update_interval:
                    prev_time = time()
                    if not infile.closed:
                        pos = infile_pos()
                    else:
                        # assuming a closed infile means the entire infile has
                        # been processed.
                        pos = infile_size
                    frac = float(pos) / infile_size
                    self._listener.notify(("PROGRESSBAR", "Alignments", frac))

                if output_alignments:
                    out.write("\t".join(map(str, v_rec)) + "\n" + \
                            "\t".join(map(str, j_rec)) + "\n")

                clone = build_clone(self._ref, v_rec, j_rec,
                        self._clone_classname)

                if clone is None:
                    continue

                seqlen = len(clone.seq)
                if seqlen < self._min_seqlen:
                    continue

                # Count base qualities in the clone (which is possible because
                # at this point the clone is based on a single read)
                lenfam_Q_counts = Q_counts.setdefault(seqlen, [0] * 42)
                for i in xrange(clone.v.end, clone.j.start):
                    lenfam_Q_counts[clone.qual[i]] += 1
                
                cs.add(clone, merge = True)

                v_allele = self._ref[v_rec.RNAME]
                j_allele = self._ref[j_rec.RNAME]
                # Count errors in the alignments
                for (rec, r_roi_start, r_roi_end) in \
                        ((v_rec, v_allele.refpos + v_refpos_offset, 0),
                        (j_rec, 0, j_allele.refpos + j_refpos_offset)):
                    allele = self._ref[rec.RNAME]
                    lenfam_alnstats = alnstats[allele.region].setdefault(
                            seqlen, {
                        "n"     : 0,
                        "mm"    : 0,
                        "ins"   : 0,
                        "dels"  : 0,
                        "Q_mm"  : [0] * 42,
                        "Q_n"   : [0] * 42})
                    n, mm, ins, dels, r_roi_as, r_roi_ae = get_error_stats(rec,
                            allele.seq,
                            lenfam_alnstats["Q_mm"], lenfam_alnstats["Q_n"],
                            r_roi_start, r_roi_end)
                    lenfam_alnstats["n"] += n
                    lenfam_alnstats["mm"] += mm
                    lenfam_alnstats["ins"] += ins
                    lenfam_alnstats["dels"] += dels
        finally:
            if output_alignments:
                out.close()
        self._listener.notify(("PROGRESSBAR", "Alignments", "end"))

        if len(cs) == 0:
            msg = "No clones found in alignments. \
Was the correct germline reference used?"
            logger.error(msg)
            raise Exception(msg)

        if not self._alignment_stats_fn is None and \
                self._alignment_stats_fn != "":
            logger.info("Writing alignment stats to \"%s\""%
                    self._alignment_stats_fn)
            with zopen(self._alignment_stats_fn, 'w') as out:
                out.write("seqlen,region,n,mm,ins,dels\n")
                for region in alnstats:
                    for seqlen, lenfam_alnstats in \
                            alnstats[region].iteritems():
                        out.write(",".join(map(str,[
                            seqlen, region,
                            lenfam_alnstats["n"],
                            lenfam_alnstats["mm"],
                            lenfam_alnstats["ins"],
                            lenfam_alnstats["dels"]])) + "\n")

        self._save_cloneset(cs, "r")

        # Sum all the counts in the V and J regions separately, and calculate
        # average error rates
        tot_err = {"V":{}, "J":{}}
        for region in ("V", "J"):
            region_stats = alnstats[region]
            x = tot_err[region]
            x["n"] = sum([y["n"] for y in region_stats.itervalues()])
            x["mm"] = sum([y["mm"] for y in region_stats.itervalues()])
            x["ins"] = sum([y["ins"] for y in region_stats.itervalues()])
            x["dels"]= sum([y["dels"] for y in region_stats.itervalues()])

            n = x["n"]
            if n > 0:
                x["mmr"] = float(x["mm"]) / n 
                x["insr"] = float(x["ins"]) / n
                x["delsr"] = float(x["dels"]) / n
            else:
                x["mmr"] = 0
                x["insr"] = 0
                x["delsr"] = 0
        global_mmr = max(tot_err["V"]["mmr"], tot_err["J"]["mmr"])
        global_insr = max(tot_err["V"]["insr"], tot_err["J"]["insr"])
        global_delsr = max(tot_err["V"]["delsr"], tot_err["J"]["delsr"])
        logger.info("global error rates: mmr: %(global_mmr)s, \
insr: %(global_insr)s, delsr: %(global_delsr)s"%locals())

        # Calculate observed error rates for Phred scores
        Q_mm_stats = {"V":{}, "J":{}}
        for region, region_stats in alnstats.iteritems():
            Q_mm = Q_mm_stats[region].setdefault("Q_mm", [0] * 42)
            Q_n = Q_mm_stats[region].setdefault("Q_n", [0] * 42)
            for lenfam_alnstats in region_stats.itervalues():
                for i in xrange(42):
                    Q_mm[i] += lenfam_alnstats["Q_mm"][i]
                    Q_n[i] += lenfam_alnstats["Q_n"][i]

        if not self._Q_mm_stats_fn is None and self._Q_mm_stats_fn != "":
            with zopen(self._Q_mm_stats_fn, 'w') as out:
                out.write("region,Q,n,mm\n")
                for region in Q_mm_stats:
                    for Q,(mm, n) in enumerate(izip(Q_mm_stats[region]["Q_mm"],
                            Q_mm_stats[region]["Q_n"])):
                        out.write("%s,%s,%s,%s\n"%(region, Q, n, mm))

        # Calculate ratio between base quality score assigned by the sequencer
        # and observed base quality (based on alignments with germline
        # reference).
        sum_ratios = 0
        n_ratios = 0
        for region in Q_mm_stats:
            Q_mm = Q_mm_stats[region]["Q_mm"]
            Q_n = Q_mm_stats[region]["Q_n"]
            for q in xrange(42):
                mm = Q_mm[q]
                n = Q_n[q]
                if mm > 0 and n > 0:
                    q_obs = p2q(float(mm) / n)
                    if q_obs > 0:
                        sum_ratios += (q / q_obs) * n
                        n_ratios += n
        if n_ratios > 0:
            alpha = float(sum_ratios) / n_ratios
        else:
            logger.warning('No instances found of a Phred score associated ' +\
                    'with mismatches.')
            alpha = 1.0

        logger.info("Ratio between base quality and observed quality: %s"%
                alpha)

        if not self._Q_mm_stats_plot_fn is None and \
                self._Q_mm_stats_plot_fn != "":
            plot_Q_mm_stats(Q_mm_stats, self._Q_mm_stats_plot_fn)

        # Get median quality score
        Q_n = [0] * 42 # count number of bases for every Q score
        for lenfam_Q_counts in Q_counts.itervalues():
            for q, count in enumerate(lenfam_Q_counts):
                Q_n[q] += count
        i = ((sum(Q_n) + 1) // 2) - 1 # index of median element in Q_n
        j = 0
        for max_Q, count in enumerate(Q_n):
            j += count
            if j > i:
                break
        logger.info("max_Q = %s"%max_Q)

        pool = ConnectedConsumerPool(n_consumers = self._n_threads)
        by_seqlen = lambda clone:len(clone.seq)
        confidence = self._confidence
        for seqlen, clones in groupby(sorted(cs, key = by_seqlen), by_seqlen):
            if self.stopped():
                logger.warning("Pipeline stopped")
                return
            cs2 = CloneSet(clones)
            # Calculate expected number of errors based on Q scores
            lenfam_Q_counts = Q_counts[seqlen]

            # get total number of bases between V and J region
            n_o = sum(lenfam_Q_counts)
            mm_o = 0
            for q, count in enumerate(lenfam_Q_counts):
                q /= alpha
                mm_o += q2p(q) * count

            mm_v = alnstats["V"][seqlen]["mm"]
            n_v = alnstats["V"][seqlen]["n"]

            mm_j = alnstats["J"][seqlen]["mm"]
            n_j = alnstats["J"][seqlen]["n"]

            mm_tot = mm_v + mm_o + mm_j
            n_tot = n_v + n_o + n_j
            logger.info("Mismatch stats for seqlen %s: mm_v (%s, %s),\
mm_o (%s, %s), mm_j (%s, %s), mm_tot (%s, %s)"%(seqlen,
                mm_v, float(mm_v)/n_v if n_v > 0 else 0,
                mm_o, float(mm_o)/n_o if n_o > 0 else 0,
                mm_j, float(mm_j)/n_j if n_j > 0 else 0,
                mm_tot, float(mm_tot)/n_tot if n_tot > 0 else 0))
            local_mmr = float(mm_tot) / n_tot
            mmr = max(local_mmr, global_mmr)
            logger.info("Adding task: seqlen: %(seqlen)s, mismatch_rate: \
%(mmr)s, confidence: %(confidence)s, max_Q: %(max_Q)s"%locals())
            pool.add_task(run_ec_on_bin, (cs2, mmr, confidence, max_Q))
    
        self._listener.notify("Running QMerge and IMerge on bins.")
        self.run_pool(pool, desc = 'QMerge, IMerge')
        results = pool.results
        cloneset = CloneSet(chain.from_iterable([x[0] for x in results]))
        self._save_cloneset(cloneset, "rqi")

        self._listener.notify("Running LMerge")
        cloneset = run_lmerge(cloneset, global_mmr, global_insr, global_delsr,
                confidence)
        self._save_cloneset(cloneset, "rqil")

        pool = ConnectedConsumerPool(n_consumers = self._n_threads)
        for seqlen, clones in groupby(sorted(cloneset, key = by_seqlen),
                by_seqlen):
            cs2 = CloneSet(clones)
            pool.add_task(wrapper_run_nmerge_on_bin,
                    args = (cs2,))
        self._listener.notify("Running NMerge on bins.")
        self.run_pool(pool, desc = 'NMerge')
        results = pool.results
        cloneset = CloneSet(chain.from_iterable(results))
        self._save_cloneset(cloneset, "rqiln")

        ########################
        # Write clones to file #
        ########################
        self._listener.notify("Writing clones")
        sequence_id = 0
        with open(self._output_fn, 'w') as res_ok:
            with open(self._output_not_ok_fn, 'w') as res_not_ok:
                header = '\t'.join(\
                        clone2AIRRDict(clone = None, ref = None).keys()) + '\n'
                res_ok.write(header)
                res_not_ok.write(header)

                n_discarded = 0
                for clone in sorted(cloneset,
                        key = lambda clone:(-clone.count, clone.seq)):
                    record = clone2AIRRDict(clone = clone, ref = self._ref)
                    min_phred = int(record['junction_minimum_quality_score'])
                    if min_phred < self._min_phred_threshold \
                            or record['stop_codon'] == 'T' or \
                            record['vj_in_frame'] == 'F':
                        n_discarded += 1
                        out = res_not_ok
                    else:
                        out = res_ok
                    sequence_id += 1
                    record['sequence_id'] = str(sequence_id)
                    out.write('\t'.join([v for k, v in record.iteritems()]) +\
                            '\n')
        self._listener.notify("Discarded %s clones"%n_discarded)

    def stopped(self):
        return self._stop.is_set()

    def stop(self):
        self._stop.set()
