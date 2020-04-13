#!/usr/bin/env python
#
# Description: performs demultiplexing and UMI (if any) extraction.

import argparse
import os.path
import gzip
from sys import stdout
from time import time

from barcode import Adapter, Barcode, revcomp
from fileio import FastqFormat, BarcodeFormat, zopen, FastqRecord, filesize
import terminal as term

def add_parser_arguments(parser):
    parser.add_argument('-f', '--forward',
            help = 'fastq file with forward reads')
    parser.add_argument('-r', '--reverse',
            help = 'fastq file with reverse reads')
    parser.add_argument('-p', '--paired',
            help = 'fastq file with reads containing merged read pairs')
    parser.add_argument('-b', '--barcodes', required = True,
            help = 'file containing sample barcodes and UMIs')
    parser.add_argument('-m', '--max_mm', type = int, default = 2,
            help = 'Number of mismatches allowed when determining length of \
region that is matching between a pair of sequences.')
    parser.add_argument('-rc', '--reverse_complement', action = 'store_true',
            help = 'Also look for (master) barcodes on reverse complement')
    return parser

def get_best_match(matches, matches_rc):
    is_rc = False
    best_match = None

    if matches and len(matches) > 0:
        best_match = matches[0]

    if matches_rc and len(matches_rc) > 0:
        if best_match is None or matches_rc[0][0] < best_match[0]:
            best_match = matches_rc[0]
            is_rc = True
    return best_match, is_rc

def checkout(fq1_fn, fq2_fn, adapters, max_mm, search_rc, paired = False):
    assert not fq1_fn is None
    assert not (paired and not fq2_fn is None)

    print 'Handling file(s): %s'%''.join([fq1_fn,
        '' if fq2_fn is None else ', %s'%fq2_fn])

    fq1_file = zopen(fq1_fn, 'r')
    if isinstance(fq1_file, gzip.GzipFile):
        fq1_filesize = os.path.getsize(fq1_file.name)
        fq1_filepos = fq1_file.fileobj.tell
    else:
        fq1_filesize = filesize(fq1_file)
        fq1_filepos = fq1_file.tell

    fq1 = FastqFormat.records_in(fq1_file, encoding = None)
    if not fq2_fn is None:
        fq2 = FastqFormat.records_in(zopen(fq2_fn, 'r'), encoding = None)
    else:
        fq2 = None

    outfiles = {}
    for (sample_id, master, slave) in adapters:
            outfiles[sample_id] = {
                    "out1" : (open("%s_R1.fastq"%sample_id, 'w'), 'R1') \
                            if not paired else \
                            (open("%s_R12.fastq"%sample_id, 'w'), 'R12'),
                    "out2" : (None, None) if fq2 is None else \
                            (open("%s_R2.fastq"%sample_id, 'w'), 'R2')}

    n_accepted = 0
    prev_time = time()
    for i, r1 in enumerate(fq1):
        if fq2:
            r2 = next(fq2)
            assert(r1.id == r2.id)
        else:
            r2 = None

        # Demultiplex
        best_match = None
        for (sample_id, master, slave) in adapters:
            matches, matches_rc = master.locate_in(r1.seq, max_mm, search_rc)

            master_match, is_rc = get_best_match(matches, matches_rc)
            
            # look for master on the mate
            if (not master_match or master_match[0] < len(master.seq)) and \
                    fq2:
                # master not found or no full length match found
                matches2, matches_rc2 = master.locate_in(r2.seq, max_mm,
                        search_rc)

                master_match2, is_rc2 = get_best_match(matches2, matches_rc2)

                if not master_match2 and not master_match:
                    # master not found on r1 nor on r2
                    continue

                if not master_match or (master_match2 and \
                        master_match2[0] < master_match[0]):
                    master_match = master_match2
                    is_rc = is_rc2
                    # apparently strands are swapped
                    r1, r2 = r2, r1

            if master_match is None:
                continue

            if is_rc:
                master_match = list(master_match)
                master_match[1] = \
                        len(r1.seq) - (master_match[1] + len(master.seq))
                master_match = tuple(master_match)
                r1 = FastqRecord(id = r1.id + " rc", desc = r1.desc,
                        seq = revcomp(r1.seq), qual_str = r1.qual_str[::-1])
                if fq2:
                    r2 = FastqRecord(id = r2.id + " rc", desc = r2.desc,
                            seq = revcomp(r2.seq),
                            qual_str = r2.qual_str[::-1])

            # Master adapter has been found, retrieve its UMI (if it has one)
            master_umi = ("", "")
            if master.has_UMI(): # get umi
                master_umi = master.get_UMI(r1.seq, r1.qual_str,
                        master_match[1])
                if master.UMI_length != len(master_umi[0]):
                    # failed to retrieve UMI from master adapter
                    continue

            # Look for slave adapter
            slave_match = None
            
            slave_umi = ("", "")
            if slave: # has slave adapter
                if paired:
                    r = r1
                else:
                    r = r2

                slave_matches, slave_matches_rc = slave.locate_in(r.seq,
                        max_mm, search_rc = False)
                slave_match = get_best_match(slave_matches, slave_matches_rc)

                if slave.has_UMI(): # get umi
                    slave_umi = slave.get_UMI(r.seq, r.qual_str,
                            slave_match[1])
                    if slave.UMI_length != len(slave_umi[0]):
                        continue

            if not best_match or best_match[0][0] > master_match[0] or \
               (best_match[0][0] == master_match[0] and \
                slave and slave_match[0] and \
                (not best_match[1] or \
                 not best_match[1][0] or \
                 best_match[1][0] > slave_match[0])):
                umi = [x + y for x,y in zip(master_umi, slave_umi)]
                best_match = (master_match, slave_match, sample_id, umi)

        if best_match:
            master_match, slave_match, sample_id, umi = best_match
            for (r, (out, typename)) in ((r1, outfiles[sample_id]["out1"]),
                    (r2, outfiles[sample_id]["out2"])):
                if not out:
                    continue
                out.write("@%s UMI:%s:%s:%s\n"%(r.id, typename, umi[0],
                    umi[1]))
                out.write("%s\n+\n%s\n"%(r.seq, r.qual_str))
            n_accepted += 1
            
        frac = float(fq1_filepos()) / fq1_filesize
        if time() - prev_time > .5 or frac == 1.0:
            prev_time = time()
            stdout.write(term.EL(2) + term.CHA(0) + \
                    "Processed %s records (%.2f%%), accepted %s (%.2f%%)"%(i + 1,
                        frac*100, n_accepted,
                        (100*float(n_accepted)/(i+1))))
            stdout.flush()
            
    stdout.write('\n')

def prog_checkout(args):
    search_rc = args.reverse_complement
    barcodes_fn = args.barcodes

    adapters = list(BarcodeFormat.records_in(open(barcodes_fn, 'r')))
    
    max_mm = args.max_mm

    if args.forward is None and args.reverse is None and args.paired is None:
        raise Exception('Require at least 1 fastq input file')

    fq1_fn = args.forward
    fq2_fn = args.reverse
    fq12_fn = args.paired

    if not fq1_fn is None:
        checkout(fq1_fn, fq2_fn, adapters, max_mm, search_rc)
    elif not fq2_fn is None:
        checkout(fq2_fn, fq1_fn, adapters, max_mm, search_rc)

    if not fq12_fn is None:
        checkout(fq12_fn, None, adapters, max_mm, search_rc, paired = True)
            
def main():
    parser = argparse.ArgumentParser()
    add_parser_arguments(parser)
    args = parser.parse_args()
    prog_checkout(args)

if __name__=='__main__':
    main()
