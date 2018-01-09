#!/usr/bin/env python
# Description: Performs unique molecular identifier (UMI) based error
# correction using simple consensus sequencing.

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

import os
import gzip
from sys import stdout, stdin
from time import time

from seq import ConsensusQSequence, get_offset
from fileio import filesize
import terminal as term

def run_umi_group_ec(records, k, moff, min_score4offset, min_score4merge,
        discarded, callback = None):
    """Run unique molecular identifier (UMI) based error correction using
    simple consensus sequencing. Returns tuple containing dictionary with UMI
    name as key and consensus sequences (ConsensusQSequence object) as value.

    :records: iterable providing QSequence or QSequence-like objects. Objects
    should contain "name", "seq", and "qual" attributes, providing sequence
    name, sequence, and a list/vector of quality scores. The "name" attribute
    should have a field like: "UMI:XXX:YYY:ZZZ" where XXX is the is an
    additional qualifier to separate out forward, reverse, and merged
    read pairs, YYY is the UMI, and ZZZ is the per base corresponding ascii
    encoded Phred(+33) score. 
    :k: window size parameter of get_offset method.
    :moff: maximum offset for the window of get_offset method.
    :min_score4offset: minimum score required to accept non-zero offset in an
    alignment between two sequences. Score is calculated according to the
    get_offset function of the seq module.
    :min_score4merge: minimum score required for a sequence to be merged with
    an existing consensus sequence having the same UMI. Otherwise another
    consensus sequence is started.
    :discarded: object that can receive QSequence-like objects retrieved from
    the records parameter, using the "+" operator (e.g. discarded += r, where
    r is a QSequence-like object).
    :callback: function to call after processing a record. Should take two
    parameters, "rec_nr" and "n_discarded", the number of discarded
    records and record number (0-based), respectively.
    """
    umi_groups = {}
    n_discarded = 0
    for rec_nr, r in enumerate(records):
        # Retrieve UMI from the name of the record
        fields = r.name.split("UMI:")
        assert len(fields) > 0
        fields = fields[1].split(":")
        assert len(fields) > 0
        qualifier, umi, qual_str = fields[:3]
        qual_str = ':'.join(fields[2:])
        assert len(qual_str) == len(umi)

        if not qualifier in ['R1', 'R2', 'R12']:
            msg = 'Format error: unknown qualifier (\'%s\') in UMI'%qualifier
            logger.error(msg)
            raise Exception(msg)
        
        name = '%s:%s'%(qualifier, umi)

        if len(r.seq) < k:
            if not discarded is None:
                discarded += r
            continue

        if name in umi_groups:
            # Due to PCR and sequencing errors, different amplicons may end up
            # with the same UMI. Here we search for the most similar sequence
            # in the current UMI group.
            comparisons = [(get_offset(r.seq, other.seq, k = k,
                    moff = moff), other) for other in umi_groups[name]]

            ((offset, score), other) = sorted(comparisons,
                    key = lambda ((offset, score), other):score,
                    reverse = True)[0]

            # Check if the current sequence should be added to the consensus of
            # 'other'
            if score < min_score4merge:
                umi_groups[name].append(ConsensusQSequence(r))
            else:
                if score < min_score4offset:
                    offset = 0
                other.add(ConsensusQSequence(r), offset)
        else:
            umi_groups[name] = [ConsensusQSequence(r)]

        if not callback is None:
            callback(rec_nr, n_discarded)
    return umi_groups

def add_parser_arguments(parser):
    parser.add_argument("-i")
    parser.add_argument("-o", required = True)
    parser.add_argument("-k", type = int, default = 30,
            help = "window width for alignment between sequences")
    parser.add_argument("-moff", "--max_offset", type = int, default = 4,
            help = "maximum offset for alignment between sequences")
    parser.add_argument("-s4o", "--min_score_for_offset", type = int,
            default = 20,
            help = "minimum score for non-zero offset in alignment between \
sequences.")
    parser.add_argument("-s4m", "--min_score_for_merge", type = int,
            default = 24,
            help = "minimum score for a sequence to merge with an existing \
consensus sequence having the same UMI.")
    return parser

def _progress_indicator(rec_nr, n_discarded, handle = None):
    if not handle is None:
        if isinstance(handle, gzip.GzipFile):
            _progress_indicator.endpos = os.path.getsize(handle.name) 
            _progress_indicator.get_pos = handle.fileobj.tell
        else:
            try:
                _progress_indicator.endpos = filesize(handle)
                _progress_indicator.get_pos = handle.tell
            except IOError:
                _progress_indicator.endpos = None
                _progress_indicator.get_pos = lambda : None

        _progress_indicator.prev_time = time()
        return
    get_pos = _progress_indicator.get_pos
    endpos = _progress_indicator.endpos

    frac = None if endpos is None else float(get_pos()) / endpos
    if time() - _progress_indicator.prev_time > .5 or frac == 1.0:
        _progress_indicator.prev_time = time()
        perc_str = '?%' if frac is None else '%.2f%%'%(frac*100)
        stdout.write(term.EL(2) + term.CHA(0) + \
                'Processed %s records (%s)'%(rec_nr + 1, perc_str))
        stdout.flush()
   
def prog_umi_group_ec(args):
    import os
    from fileio import zopen, FastqFormat
    from seq import QSequence

    assert args.k > 0
    assert args.max_offset >= 0
    assert args.min_score_for_offset <= args.k
    assert args.min_score_for_merge <= args.k

    fn, fext = os.path.splitext(args.o)
    discarded = FastqFormat.records_out(
            zopen(fn + '.discarded' + fext, 'w'), None)

    input_file = stdin if args.i is None else zopen(args.i, 'r')
    _progress_indicator(0, 0, input_file)
    umi_groups = run_umi_group_ec(
            records = FastqFormat.records_in(input_file),
            k = args.k,
            moff = args.max_offset,
            min_score4offset = args.min_score_for_offset,
            min_score4merge = args.min_score_for_merge,
            discarded = discarded,
            callback = _progress_indicator)
    stdout.write('\n')

    out = FastqFormat.records_out(zopen(args.o, 'w'), None)

    print 'Writing results to file'
    for name, grouplist in umi_groups.iteritems():
        grouplist = sorted(grouplist, key = lambda cqs:(-cqs.count,
            -sum(cqs.qual)))
        cqs = grouplist[0]
        out += QSequence('UMI:%s:%s:%s'%(name,
            '%s/%s'%('1', len(grouplist)), cqs.count), cqs.seq, cqs.qual)
        for i, cqs in enumerate(grouplist[1:]):
            discarded += QSequence('UMI:%s:%s:%s'%(name, '%s/%s'%(i+2,
                len(grouplist)),cqs.count),
                    cqs.seq, cqs.qual)
        
def main():
    import argparse
    parser = argparse.ArgumentParser()
    add_parser_arguments(parser)
    args = parser.parse_args()

    prog_umi_group_ec(args)

if __name__ == "__main__":
    main()
