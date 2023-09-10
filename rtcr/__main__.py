__version__ = "0.5.2"

import logging
import logging.config
import os
from os import path
import sys
from sys import stdout
import argparse
import multiprocessing as mp
from collections import OrderedDict

import terminal as term
import config
import convert

# NOTE: logger should not be used until after logging has been initialised.
logger = logging.getLogger(__name__)
here = path.abspath(path.dirname(__file__))

def handle_uncaught_exception(exc_type, exc_value, exc_traceback):
    logger.critical('Uncaught exception', exc_info = (exc_type, exc_value,
        exc_traceback))

def init_logging():
    logging.config.fileConfig(config.config_fn,
            disable_existing_loggers = False)

def determine_phred_encoding(reads_fn, n = 100):
    with open(reads_fn,'r') as f:
        for i, line in enumerate(f):
            if i >= n*4:
                break
            if i % 4 == 3:
                for qch in line.rstrip('\n'):
                    oqch = ord(qch)
                    if oqch < 59:
                        return 33
    return 64

class Listener(object):
    def __init__(self):
        self.progressbars = OrderedDict() 

    def notify(self, msg):
        if isinstance(msg, basestring):
            print msg
        elif msg[0] == "PROGRESSBAR":
            n = len(self.progressbars)
            name, action = msg[1:]
            if action == "start":
                stdout.write(term.CLRSCR())
                self.progressbars[name] = term.ProgressBar()
            elif action == "end":
                stdout.write(term.CLRSCR())
                del self.progressbars[name]
            else:
                self.progressbars[name].frac = action
            stdout.write(term.CUP(0, 0))
            for name, progressbar in self.progressbars.iteritems():
                stdout.write("%s : %s\n"%(progressbar, name))

def str2bool(s):
    """Convert string to boolean.

    :s: string to convert
    """
    if s.lower() == "true":
        return True
    elif s.lower() == "false":
        return False
    else:
        raise ValueError("string not matching 'true' or 'false'")

def pipeline_add_parser_arguments(parser):
    parser.add_argument("--ref", help = "Immune receptor reference file")
    parser.add_argument("-i","--reads", required = True,
            help = "Fastq file with reads")
    parser.add_argument("-p","--phred_encoding", type = int,
            help = "Phred encoding (ascii code for a PHRED score of 0)")
    parser.add_argument("-t", "--threads", type = int,
            help = "Number of threads or processes")
    parser.add_argument("--species",help = "Reference species; use comma to \
select multiple species.")
    parser.add_argument("--gene", help = "Reference gene; use comma to \
select multiple genes.")
    parser.add_argument("--no_VJ_collapsing", action = "store_true",
            help = "Do not collapse clones with identical CDR3 and different \
VJ combination.")
    parser.add_argument("--debug", action = "store_true",
            help = "Set logging level to DEBUG")
    parser.add_argument("-v", "--verbose", action = "store_true",
            help = "Output log statements to stdout")
    return parser

def parse_cmdline():
    import checkout
    import umi_group_ec

    # top-level parser
    parser = argparse.ArgumentParser(prog = "rtcr")
    
    subparsers = parser.add_subparsers()

    # Add Config program
    parser_Config = subparsers.add_parser("Config",
            help = "Configure pipeline")
    parser_Config = config.add_parser_arguments(parser_Config)
    parser_Config.set_defaults(func = config.prog_config)

    # Add Checkout program
    parser_Checkout = subparsers.add_parser("Checkout",
            help = "Perform demultiplexing and UMI extraction",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser_Checkout = checkout.add_parser_arguments(parser_Checkout)
    parser_Checkout.set_defaults(func = checkout.prog_checkout)

    # Add UMI_group_ec program
    parser_umi_group_ec = subparsers.add_parser("umi_group_ec",
            help = "Perform UMI group error correction",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser_umi_group_ec = umi_group_ec.add_parser_arguments(
            parser_umi_group_ec)
    parser_umi_group_ec.set_defaults(func = umi_group_ec.prog_umi_group_ec)

    # Add convert program
    parser_Convert = subparsers.add_parser('Convert',
            help = 'Convert cloneset (.dat) files to human readable format')
    parser_Convert = convert.add_parser_arguments(parser_Convert)
    parser_Convert.set_defaults(func = convert.prog_Convert)

    # Add pipeline program
    parser_Pipeline = pipeline_add_parser_arguments(subparsers.add_parser(
        "run", help = "main rtcr pipeline"))
    parser_Pipeline.set_defaults(func = prog_pipeline)

    if len(sys.argv) == 1:
        print "RTCR version %s"%__version__
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    args.func(args)

def prog_pipeline(args):
    # [Defaults] section
    confidence = float(config.get("Defaults", "confidence"))
    min_phred_threshold = int(config.get("Defaults", "min_phred_threshold"))
    clone_classname = config.get("Defaults", "clone_classname")
    min_seqlen = int(config.get("Defaults", "min_seqlen"))
    species = config.get("Defaults", "species")
    gene = config.get("Defaults", "gene")
    update_interval = float(config.get("Defaults", "update_interval"))
    ref_fn = path.join(here, config.get("Defaults", "ref_fn"))
    alignments_fn = config.get("Defaults", "alignments_fn")
    alignment_stats_fn = config.get("Defaults", "alignment_stats_fn")
    Q_mm_stats_fn = config.get("Defaults", "Q_mm_stats_fn")
    Q_mm_stats_plot_fn = config.get("Defaults", "Q_mm_stats_plot_fn")
    qplot_fn = config.get("Defaults", "qplot_fn")
    output_fn = config.get('Defaults', 'output_fn')
    output_not_ok_fn = config.get('Defaults', 'output_not_ok_fn')

    # [Aligner] section
    location = config.get("Aligner", "location")
    cmd_build_index = config.get("Aligner", "cmd_build_index")
    args_build_index = config.get("Aligner", "args_build_index")
    cmd_align = config.get("Aligner", "cmd_align")
    args_align_base = config.get("Aligner", "args_align_base")
    args_align_v = args_align_base + " " + \
            config.get("Aligner", "args_align_v")
    args_align_j = args_align_base + " " + \
            config.get("Aligner", "args_align_j")

    if args.no_VJ_collapsing:
        clone_classname = "AnnotatedCloneDistinctAllele"

    ref_fn = path.realpath(ref_fn) \
            if args.ref is None else path.realpath(args.ref)
    reads_fn = path.realpath(args.reads)
    phred_encoding = determine_phred_encoding(reads_fn) \
            if args.phred_encoding is None else args.phred_encoding
    assert phred_encoding in (33, 64)
    species =  species if args.species is None else args.species
    species = species.split(",")
    gene =  gene if args.gene is None else args.gene
    gene = gene.split(",")
    n_threads = mp.cpu_count() if args.threads is None else args.threads

    # Prepare aligner commands and check existence of aligner
    cmd_build_index = path.realpath(path.join(location, cmd_build_index))
    cmd_align = path.realpath(path.join(location, cmd_align))

    if not path.isfile(cmd_build_index):
        raise ValueError(\
                "Executable to build an index (\"%s\") does not exist.\n\
Please use \"rtcr Config\" to see if the Aligner is properly configured"%
cmd_build_index)

    if not path.isfile(cmd_align):
        raise ValueError(\
                "Executable to align sequences (\"%s\") does not exist.\n\
Please use \"rtcr Config\" to see if the Aligner is properly configured "%
cmd_align)

    init_logging()
    if args.debug:
        logging.root.setLevel(logging.DEBUG)
        logging.root.debug("log level set to DEBUG")
    if args.verbose:
        logging.root.addHandler(logging.StreamHandler(stdout,))
        logging.root.info("Writing log statements to stdout")
    # Note, delaying import of modules that have a logger until after logging
    # has been initialised.
    from fileio import read_reference, zopen
    from pipeline import Pipeline

    ref = read_reference(ref_fn).get_slice(species = species, gene = gene)

    for s in species:
        if not s in ref.species:
            logger.error("species \"%s\" does not exist in reference"%s)
            sys.exit(1)
    for g in gene:
        if not g in ref.genes:
            logger.error("gene \"%s\" does not exist in reference"%s)
            sys.exit(1)

    version = __version__
    preamble = '\nRTCR version: %(version)s\n'%locals()
    preamble += '\n[Command line arguments]\n' + \
            '\n'.join(['%s : %s'%(i,v) for i,v in enumerate(sys.argv)])
    preamble += '\n\
[Files]\n\
Reference: %(ref_fn)s\n\
Reads: %(reads_fn)s\n\
Output: %(output_fn)s\n\
\n\
[Settings]\n\
PHRED encoding: %(phred_encoding)s\n\
Species: %(species)s\n\
Gene: %(gene)s\n\
confidence: %(confidence)s\n\
\n\
[Immune receptor reference]\n'%locals()
    for species in sorted(ref.species):
        for gene in sorted(ref.genes):
            for region in sorted(ref.regions):
                alleles = ref.get_alleles(species = species, gene = gene,
                        region = region)
                n = len(alleles)
                if n > 0:
                    preamble += "%s,%s,%s: %s alleles\n"%(
                            species, gene, region, n)
                s = ""
                for allele in alleles:
                    s += "%s, %s\n"%(allele.species, allele.name)
                logger.debug("species, allele\n" + s)
    preamble += "\n[Pipeline run]"
    logger.info(preamble)

    # Make sure exceptions are logged, even when not caught
    sys.excepthook = handle_uncaught_exception

    pipeline = Pipeline(
            ref = ref,
            reads = zopen(reads_fn, 'r'),
            phred_encoding = phred_encoding,
            cmd_build_index = cmd_build_index,
            args_build_index = args_build_index,
            cmd_align = cmd_align,
            args_align_v = args_align_v,
            args_align_j = args_align_j,
            alignments_fn = alignments_fn,
            alignment_stats_fn = alignment_stats_fn,
            Q_mm_stats_fn = Q_mm_stats_fn,
            Q_mm_stats_plot_fn = Q_mm_stats_plot_fn,
            output_fn = output_fn,
            output_not_ok_fn = output_not_ok_fn,
            clone_classname = clone_classname,
            confidence = confidence,
            min_seqlen = min_seqlen,
            min_phred_threshold = min_phred_threshold,
            n_threads = n_threads,
            update_interval = update_interval,
            listener = Listener())
    pipeline.daemon = True
    pipeline.name = 'Pipeline'
    try:
        pipeline.start()
        while pipeline.is_alive():
            pipeline.join(1)
    except KeyboardInterrupt:
        logger.error('Caught keyboard interrupt. Shutting down.')
        pipeline.stop()
        pipeline.join(1)

def main():
    parse_cmdline()
