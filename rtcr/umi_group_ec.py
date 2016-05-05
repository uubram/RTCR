#!/usr/bin/env python
# Description: Performs unique molecular identifier (UMI) based error
# correction using simple consensus sequencing.

import logging
logger = logging.getLogger(__name__)

from seq import ConsensusQSequence, get_offset

def run_umi_group_ec(records, k, moff, min_score, discarded = None):
    """Run unique molecular identifier (UMI) based error correction using
    simple consensus sequencing. Returns tuple containing dictionary with UMI
    as key and consensus sequences (ConsensusQSequence object) as value.

    :records: iterable providing QSequence or QSequence-like objects. Objects
    should contain "name", "seq", and "qual" attributes, providing sequence
    name, sequence, and a list/vector of quality scores. The "name" attribute
    should have a field like: "UMI:XXX:XXX" where the first XXX is the UMI and
    the second XXX is the per base corresponding ascii encoded Phred(+33)
    score.
    :k: window size parameter of get_offset method.
    :moff: maximum offset for the window of get_offset method.
    :min_score: minimum score required to accept non-zero offset in an
    alignment between two sequences. Score is calculated according to the
    get_offset function of the seq module.
    :discarded: object that can receive QSequence-like objects retrieved from
    the records parameter, using the "+" operator (e.g. discarded += r, where
    r is a QSequence-like object).
    """
    umi_groups = {}
    for r in records:
        # Retrieve UMI from the name of the record
        fields = r.name.split("UMI:")
        assert len(fields) > 0
        fields = fields[1].split(":")
        assert len(fields) > 0
        umi, qual_str = fields[:2]

        if len(r.seq) < k:
            if not discarded is None:
                discarded += r
            continue

        if umi in umi_groups:
            offset, score = get_offset(r.seq, umi_groups[umi].seq, k = k,
                    moff = moff)
            if score < min_score:
                offset = 0
            umi_groups[umi].add(ConsensusQSequence(r), offset)
        else:
            umi_groups[umi] = ConsensusQSequence(r)
    return umi_groups

def add_parser_arguments(parser):
    parser.add_argument("-i", required = True)
    parser.add_argument("-o", required = True)
    parser.add_argument("-k", type = int, default = 30,
            help = "window width for alignment between sequences")
    parser.add_argument("-moff", "--max_offset", type = int, default = 4,
            help = "maximum offset for alignment between sequences")
    parser.add_argument("-s", "--min_score", type = int, default = 20,
            help = "minimum score for non-zero offset in alignment between \
sequences.")
    return parser

def prog_umi_group_ec(args):
    import os
    from fileio import zopen, FastqFormat
    from seq import QSequence

    assert args.k > 0
    assert args.max_offset >= 0

    fn, fext = os.path.splitext(args.o)
    with zopen(fn + ".discarded" + fext, 'w') as discard_file:
        umi_groups = run_umi_group_ec(FastqFormat.records_in(
            zopen(args.i, 'r')), args.k, args.max_offset, args.min_score,
            FastqFormat.records_out(discard_file, None))

    FastqFormat.records_out(zopen(args.o, 'w'),
            (QSequence("UMI:%s:%s"%(umi, cqs.count),
                cqs.seq, cqs.qual) for umi, cqs in umi_groups.iteritems()))

def main():
    import argparse
    parser = argparse.ArgumentParser()
    add_parser_arguments(parser)
    args = parser.parse_args()

    prog_umi_group_ec(args)

if __name__ == "__main__":
    main()
