#!/usr/bin/env python

#
# Copyright (C) 2006 Victoria University of Wellington
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
Starts an FSindex search server.

SYNOPSIS:

    FSsearchc.py port indexfile

DESCRIPTION:

This is non-daemonic, platform independent version of FSsearchd.
Logs everything to the standard output.

ARGUMENTS:

port:      Port that the server will bind.
indexfile: Path to the index to load.

If the inxdexfile ends with .cfg it is treated as a configuration
file for slaves. Each line should be in the following format:
host port workpath indexfile [pythonpath [binpath]]

Note that only host and port parameters are used since the master
does not start the slaves automatically.
"""

import sys, signal
from ShortFrags.Expt import SearchServer

__version__ = "$Rev: 1.1 $"
__date__ = "$Date: 2006-03-02 15:48:06 -0500 (Thu, 02 Mar 2006) $"
__author__ = "Aleksandar Stojmirovic"

SrchS = None

def signal_handler(signum, frame):
    global SrchS
    SrchS.terminate_flag = SearchServer.SIGNAL

if __name__=='__main__':
    if len(sys.argv) < 3:
        print __doc__
        sys.exit(1)

    port = int(sys.argv[1])
    indexfile = sys.argv[2]

    SrchS = SearchServer.SearchServer('CONSOLE', port, indexfile, False)
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)
    SrchS.start()
