#!/usr/bin/env python
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
__date__ = "$Date: 2005/07/13 07:01:27 $"
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




