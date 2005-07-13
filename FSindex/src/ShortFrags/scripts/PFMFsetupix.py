#!/usr/bin/env python
"""
Creates FSindexes according to a configuration XML file.

SYNOPSIS:

    PFMFsetupix [config_file]
    PFMFsetupix -h|--help


DESCRIPTION:

This script creates instances of FSindex according to the alphabet partitions
given in the configuration XML file. The XML file can be specified as the 
command-line argument or entered through standard input.
"""

import getopt, sys
from ShortFrags.Setup.SetupIndex import PFMF_IndexCreator

__version__ = "$Revision: 1.1 $"
__date__ = "$Date: 2005/07/13 06:57:57 $"
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

    IL = PFMF_IndexCreator()
    IL.parse_config(fp)
    fp.close()

    IL.create_fasta_files()
    IL.create_indexes()
