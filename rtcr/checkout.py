#!/usr/bin/env python
#
# Description: performs demultiplexing and UMI (if any) extraction.

import argparse
import os.path
import gzip
from sys import stdout
from barcode import Adapter, Barcode, revcomp
from fileio import FastqFormat, BarcodeFormat, zopen, FastqRecord, filesize
import terminal as term
from time import time

def add_parser_arguments(parser):
    parser.add_argument("-i", required = True,
            help = "fastq file with reads containing UMIs")
    parser.add_argument("-i2",
            help = "second fastq file. Only required in case of paired-end \
sequencing.")
    parser.add_argument("-b", "--barcodes", required = True,
            help = "file containing sample barcodes and UMIs")
    parser.add_argument("-m", "--max_mm", type = int, default = 2,
            help = "Number of mismatches allowed when determining length of \
region that is matching between a pair of sequences.")
    parser.add_argument("-rc", "--reverse_complement", action = "store_true",
            help = "Also look for (master) barcodes on reverse complement")
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

def prog_checkout(args):
    search_rc = args.reverse_complement
    barcodes_fn = args.barcodes

    adapters = list(BarcodeFormat.records_in(open(barcodes_fn, 'r')))
    
    outfiles = {sample_id : open("%s.fastq"%sample_id,'w') \
            for (sample_id, master, slave) in adapters}

    fq1_fn = args.i
    fq2_fn = args.i2

    max_mm = args.max_mm

    fq1_file = zopen(fq1_fn, 'r')
    if isinstance(fq1_file, gzip.GzipFile):
        fq1_filesize = os.path.getsize(fq1_file.name)
        fq1_filepos = fq1_file.fileobj.tell
    else:
        fq1_filesize = filesize(fq1_file)
        fq1_filepos = fq1_file.tell

    fq1 = FastqFormat.records_in(fq1_file, encoding = None)
    if fq2_fn:
        fq2 = FastqFormat.records_in(zopen(fq2_fn, 'r'), encoding = None)
    else:
        fq2 = False

    n_accepted = 0
    prev_time = time()
    for i, r1 in enumerate(fq1):
        if time() - prev_time > .5:
            prev_time = time()
            frac = float(fq1_filepos()) / fq1_filesize
            stdout.write(term.EL(2) + term.CHA(0) + \
                    "Processed %s records (%.2f%%), accepted %s (%.2f%%)"%(i,
                        frac*100, n_accepted, 100*float(n_accepted)/i))
            stdout.flush()

        if fq2:
            r2 = next(fq2)
            assert(r1.id == r2.id)

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

            if master_match == None:
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
                slave_matches, slave_matches_rc = slave.locate_in(r2.seq,
                        max_mm, search_rc = False)
                slave_match = get_best_match(slave_matches, slave_matches_rc)

                if slave.has_UMI(): # get umi
                    slave_umi = slave.get_UMI(r2.seq, r2.qual_str,
                            slave_match[1])
                    if slave.UMI_length != len(slave_umi[0]):
                        continue

            if not best_match or best_match[0][0] > master_match[0]:
                umi = [x + y for x,y in zip(master_umi, slave_umi)]
                best_match = (master_match, sample_id, umi)

        if best_match:
            master_match, sample_id, umi = best_match
            out = outfiles[sample_id]
            out.write("@%s UMI:%s:%s\n"%(r1.id, umi[0], umi[1]))
            out.write("%s\n+\n%s\n"%(r1.seq, r1.qual_str))
            n_accepted += 1

def main():
    parser = argparse.ArgumentParser()
    add_parser_arguments(parser)
    args = parser.parse_args()
    prog_checkout(args)

if __name__=="__main__":
    main()
