#!/usr/bin/env python
"""
Sets up a BioSQL database according to a configuration XML file.

SYNOPSIS:

    PFMFsetupdb [config_file]
    PFMFsetupdb -h|--help


DESCRIPTION:

This script creates and initialises a BioSQL schema in a relational database
(only tested with PostgreSQL) and fills it with a dataset in Uniprot
(SwissProt) format. It is also possible to enter the Uniref clusters and
InterPro domains as sequence features (annotations).

All configuration options are specified in the XML file specified as the sole
argument. Alternatively, it can be entered through standard input.
"""

import getopt, sys
from ShortFrags.Setup.SetupDb import PFMF_DatabaseLoader

__version__ = "$Revision: 1.1 $"
__date__ = "$Date: 2005/07/13 06:58:23 $"
__author__ = "Aleksandar Stojmirovic"
__credits__ = "BioSQL project, Biopython project"


if __name__ == "__main__":
    options = 'h'
    long_options = ['help']
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], options,
                                   long_options)
    except getopt.GetoptError:
        # print help information and exit:
        print __doc__
        sys.exit(2)

    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit()

    if len(args):
        fp = file(args[0])
    else:
        fp = sys.stdin

    DBL = PFMF_DatabaseLoader()
    DBL.parse_config(fp)
    fp.close()
    DBL.load_database()

