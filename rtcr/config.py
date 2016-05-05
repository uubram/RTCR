#!/usr/bin/env python

import ConfigParser
import os
import sys

here = os.path.abspath(os.path.dirname(__file__))
config_fn = os.path.join(here, "config.ini")

if not os.path.isfile(config_fn):
    raise Exception("Configuration file \"%s\" does not exist"%config_fn)

# NOTE: using RawConfigParser to prevent interpolation (some options such
# as in the Aligner section, may reference local variables rather than
# other variables in the config file)
config = ConfigParser.RawConfigParser()
config.optionxform = str # case sensitive loading of config
config.read(config_fn)

def add_parser_arguments(parser):
    parser.add_argument("option", nargs = "?",
            help = "set or get configuration option. \
Format is <section>.<key>=<value>")
    return parser

def save_config():
    with open(config_fn, 'w') as f:
        f.write("# Recover T Cell Receptor (RTCR) pipeline configuration\n")
        config.write(f)

def prog_config(args):
    if args.option is None:
        config.write(sys.stdout)
    else:
        if "=" in args.option: # Assign value to option
            location, value = args.option.split("=")
            section, option = location.split(".")
            config.set(section, option, value)
            save_config()
        else: # Get option value
            section, option = args.option.split(".")
            print config.get(section, option)

def get(*args, **kwargs):
    return config.get(*args, **kwargs)
