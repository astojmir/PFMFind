#!/usr/bin/env python

#
# Copyright (C) 2004-2006 Victoria University of Wellington
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
Runs the PFMFind GUI.

SYNOPSIS:

    PFMFind [--console] [config_file]

DESCRIPTION:

This script runs the PFMFind GUI. If config_file (XML format) is
provided, it attempts to initialise the client according to
it. Otherwise, the user will have to enter configuration arguments to
a GUI form. Optionally, a Python console can be provided.
"""

import Pmw, sys, getopt
from ShortFrags.GUI.PFMFindGUI import PFMFindGUI

__version__ = "$Revision: 1.1 $"
__date__ = "$Date$"
__author__ = "Aleksandar Stojmirovic"


if __name__ == "__main__":
    config_file = None
    globals_dict = None

    # Parse Arguments
    options = ''
    long_options = ['console']
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], options,
                                   long_options)
    except getopt.GetoptError:
        # print help information and exit:
        print __doc__
        sys.exit(2)

    for o, a in opts:
        if o in ("--console"):
            globals_dict=globals()

    if len(args):
        config_file = args[0]

    # Start GUI
    root = Pmw.initialise()
    root.title("PFMFind")
    PFMF_GUI = PFMFindGUI(root, config_file, globals_dict)
    root.mainloop()




