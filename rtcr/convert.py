#!/usr/bin/env python
# Description: Converts .dat files to human readable format 

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

from sys import stdout
import cPickle as pickle

from util import clone2str, clone2AIRRDict
import config

def add_parser_arguments(parser):
    parser.add_argument('-i', required = True,
            help = 'Input cloneset (.dat) file')
    parser.add_argument('-o',
            help = 'Output tsv file')
    return parser

def prog_Convert(args):
    infile = stdin if args.i is None else open(args.i, 'r')
    outfile = stdout if args.o is None else open(args.o, 'w')

    cloneset = pickle.load(infile)

    outfile.write('\t'.join(
        clone2AIRRDict(clone = None, ref = None).keys()) + '\n')
    for clone in sorted(cloneset, key = lambda clone:clone.count,
            reverse = True):
        outfile.write('\t'.join(
            [v for k, v in \
                    clone2AIRRDict(clone = clone, ref = None).iteritems()]) \
                    + '\n')

def main():
    import argparse
    parser = argparse.ArgumentParser()
    add_parser_arguments(parser)
    args = parser.parse_args()

    prog_Convert(args)

if __name__ == '__main__':
    main()
