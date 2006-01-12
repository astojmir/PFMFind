#!/usr/bin/env python

#
# Copyright (C) 2005-2006 Victoria University of Wellington
#
# This file is part of the PFMFind module.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2,
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#


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
__date__ = "$Date$"
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
