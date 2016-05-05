#!/usr/bin/env python

import logging
logger = logging.getLogger(__name__)

from subprocess import check_call, Popen, PIPE, STDOUT
import os
import threading as th
from Queue import Queue
import shlex
import shutil

from fileio import SAMRecord, string2SAMRecord, FastaFormat
from util import TemporaryDirectory
from seq import cigar_intervals

try:
    from subprocess import DEVNULL # Py3k
except ImportError:
    import os
    DEVNULL = open(os.devnull, "wb")

def _copyfileobj(src, dest):
    with src, dest:
        shutil.copyfileobj(src, dest)

def start_aligner(dirname, ref, region, cmd_build_index,
        args_build_index, cmd_align, args_align, phred_encoding, n_threads):
    """Starts external aligner.

    :dirname: name of directory to put temporary files such as reference
    sequences
    :ref: AlleleContainer-like containing reference sequences
    :region: string, "V" or "J" region
    :cmd_build_index: command that will be run to build an index of the
    reference sequences. If empty "" or None, this command will not be run.
    :args_build_index: string, arguments that will be provided to
    cmd_build_index
    :cmd_align: string, command that will be run to start the aligner
    :args_align: string, arguments that will be provided to cmd_align
    """
    assert region == "V" or region == "J"
    ref_fn = os.path.join(dirname, "%s_ref.fa"%region)
    index_fn = os.path.join(dirname, "%s_index"%region)
    FastaFormat.records_out(open(ref_fn, 'w'),
            ref.get_alleles(region = region))
    cmd = [cmd_build_index] + shlex.split(args_build_index%locals())
    logger.info("Running command: \"%s\""%" ".join(cmd))
    check_call(cmd, stdout = DEVNULL, stderr = STDOUT)
    cmd = [cmd_align] + shlex.split(args_align%locals())
    logger.info("Running command: \"%s\""%" ".join(cmd))
    return Popen(cmd, bufsize = -1, executable = None, stdin = PIPE,
            stdout = PIPE, stderr = DEVNULL) 

def trim_and_send(v_aligner, j_aligner, results):
    """Trims V-region from reads and sends the remainder to the aligner for
    J-region alignment.

    :v_aligner: Popen object, interface to aligner process for V alignments
    :j_aligner: Popen object, interface to aligner process for J alignments
    :results: Queue, any V alignments will be turned into SAMRecord objects and
    put on this Queue.
    """
    with j_aligner.stdin as j_in:
        for line in iter(v_aligner.stdout.readline, ''):
            rec = string2SAMRecord(line)
            if rec is None:
                continue
            results.put(rec)
            if rec.FLAG & 4:
                trimmed_read = (rec.QNAME, "", "")
            else:
                qas, qae, tas, tae = cigar_intervals(rec.CIGAR, rec.POS - 1)
                trimmed_read = (rec.QNAME, rec.SEQ[qae:],
                        rec.QUAL[qae:])
            j_in.write("@%s\n%s\n+\n%s\n"%trimmed_read)

def collect_vj_output(v_records, j_aligner, results):
    """Collects V and J alignments and stores the corresponding SAMRecords in
    2-tuples in a Queue.

    Since the V and J aligner processes are running simultaneously, the results
    of the alignments are likely to arrive out of order. This method collects
    the results as they arrive and as soon as V and J alignments to the same
    read are found, their results are put onto the results queue.

    :v_records: Queue, providing SAMRecord objects corresponding to V
    alignments.
    :j_aligner: Popen object, interface to aligner that performs J alignments.
    :results: Queue, will receive 2-tuples (v_rec, j_rec), where v_rec is the
    SAMRecord corresponding to the V alignment and j_rec is the SAMRecord to
    the J alignment.
    """
    tmp_v_recs = {} # temporarily unmatched V SAM records
    tmp_j_recs = {} # temporarily unmatched J SAM records
    for line in iter(j_aligner.stdout.readline, ''):
        j_rec = string2SAMRecord(line)
        if j_rec is None:
            continue
        tmp_j_recs[j_rec.QNAME] = j_rec

        v_rec = v_records.get()
        tmp_v_recs[v_rec.QNAME] = v_rec

        while len(tmp_v_recs) > 0:
            key = next(tmp_v_recs.iterkeys())
            if key in tmp_j_recs:
                results.put( (tmp_v_recs.pop(key), tmp_j_recs.pop(key)) )
            else:
                break

def get_vj_alignments(ref, reads, cmd_build_index, args_build_index,
        cmd_align, args_align_v, args_align_j, phred_encoding, n_threads):
    """Align V and J germline reference sequences to reads (in FastQ format).
    Yields (v_rec, j_rec) as 2-tuples, where v_rec is a SAMRecord object
    containing an alignment of a V sequence to a read identified by
    v_rec.QNAME. j_rec is similar, but contains an alignment of a J sequence to
    the trimmed version of the same read.

    :ref: AlleleContainer containing both V and J germline reference alleles.
    :reads: handle to FastQ file 
    :cmd_build_index: command that will be run to build an index of the
    reference sequences. If empty "" or None, this command will not be run.
    :args_build_index: string, arguments that will be provided to
    cmd_build_index
    :cmd_align: string, command that will be run to start the aligner
    :args_align_v: string, arguments that will be provided to cmd_align to
    align V sequences to the reads.
    :args_align_j: string, arguments that will be provided to cmd_align to
    align J sequences to the reads.
    :phred_encoding: string passed to the aligner to tell it which encoding to
    use. (e.g. for bowtie2 "33" or "64" is used).
    :n_threads: maximum number of threads/processes the aligner is allowed to
    use.
    """
    with TemporaryDirectory() as dirname:
        logger.info("Created temporary directory: \"%s\""%dirname)
        # e.g. if n_threads is 7, 4 will be given to v aligner and 3 to the
        # j aligner. This is because the v aligner is first and has to do more
        # work.
        n_threads_v = int(round(n_threads / float(2)))
        n_threads_j = n_threads // 2
        v_aligner = start_aligner(dirname, ref, "V", cmd_build_index,
            args_build_index, cmd_align, args_align_v, phred_encoding,
            n_threads_v)
        j_aligner = start_aligner(dirname, ref, "J", cmd_build_index,
            args_build_index, cmd_align, args_align_j, phred_encoding,
            n_threads_j)

        v_aligner_input_thread = th.Thread(name = "input to v_aligner",
                target = _copyfileobj,
                args = (reads, v_aligner.stdin))
        v_aligner_input_thread.start()

        # Trim V-region from reads and send trimmed reads to aligner, for J
        # segment alignment, in separate thread
        v_records = Queue()
        j_aligner_input_thread = th.Thread(name = "input to j_aligner",
                target = trim_and_send,
                args = (v_aligner, j_aligner, v_records))
        j_aligner_input_thread.start()

        vj_records = Queue()
        vj_output_thread = th.Thread(name = "collect vj output",
                target = collect_vj_output,
                args = (v_records, j_aligner, vj_records))
        vj_output_thread.start()

        while vj_output_thread.is_alive() or not vj_records.empty():
            while not vj_records.empty():
                yield vj_records.get()

def main():
    from rtcr_io import read_reference, zopen
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--ref", required = True,
            help = "Immune receptor reference file")
    parser.add_argument("--species",
            help = "Which species should be used as reference")
    parser.add_argument("--gene",
            help = "Which gene should be used as reference")
    parser.add_argument("-i","--reads", required = True,
            help = "Fastq file with reads")
    parser.add_argument("-p","--phred_encoding", type=int, required = True,
            help = "Phred encoding (ascii code for a PHRED score of 0)")
    parser.add_argument("--threads", type = int, default = 1,
            help = "Number of threads the aligner should use")
    args = parser.parse_args()

    cmd_build_index = "bowtie2-build"
    args_build_index = r"%(ref_fn)s %(index_fn)s"

    cmd_align = "bowtie2"
    args_align_v = r"-D 20 -R 3 -N 0 -i S,1,0.50 -L 8 " + \
        r"--local -x %(index_fn)s - --phred%(phred_encoding)s " + \
        r"--threads %(n_threads)s"
    args_align_j = r"-D 20 -R 3 -N 0 -i S,1,0.50 -L 8 " + \
        r"--local -x %(index_fn)s - --phred%(phred_encoding)s " + \
        r"--threads %(n_threads)s"

    ref_fn = os.path.realpath(args.ref)
    species = args.species
    gene = args.gene
    reads_fn = os.path.realpath(args.reads)
    phred_encoding = args.phred_encoding
    n_threads = args.threads

    ref = read_reference(ref_fn).get_slice(species = species,gene = gene)
    for v_rec, j_rec in get_vj_alignments(ref,
            zopen(reads_fn, "rb"), cmd_build_index,
            args_build_index, cmd_align, args_align_v, args_align_j,
            phred_encoding, n_threads):
        print "v_rec = %s"%map(str, v_rec)
        print "j_rec = %s"%map(str, j_rec)
    
if __name__ == "__main__":
    main()
