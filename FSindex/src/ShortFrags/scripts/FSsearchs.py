#!d:\Python24\python.exe

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

import sys, os.path
import win32serviceutil, win32service
from pfmfind.Expt import SearchServer

__version__ = "$Rev: 1.1 $"
__date__ = "$Date: 2006-03-02 15:48:06 -0500 (Thu, 02 Mar 2006) $"
__author__ = "Aleksandar Stojmirovic"

def start_usage():
    print \
"""
Arguments for 'start' command: workpath serverid port indexfile
workpath:  the directory where the log files are to be written.
           Must be a valid path.
serverid:  arbitrary identifier that is appended to the name of the log file.
port:      port that the server will bind.
indexfile: Path to the index to load.

If the inxdexfile ends with .cfg it is treated as a configuration file for slaves.
Each line should be in the following format:
host port workpath indexfile [pythonpath [binpath]]

Note that only host and port parameters are used by the Windows version at this
moment since the master does not start the slaves automatically.
"""    
    
def check_start_command_line():
    """
    If the 'start' command is given, it makes sure that
    1. There are at least 4 arguments, and
    2. The first argument (workdir) is a valid path.

    The reason for this is that we want to be sure that the
    messages will be written to a log so that the service
    does not die silently.
    """


    # Scan the command line in the same way as win32serviceutil.HandleCommandLine does
    import getopt
    long_opts = ["password=","username=","startup=","perfmonini=", "perfmondll=", "interactive", "wait="]
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", long_opts)
    except getopt.error, details:
        print details
        win32serviceutil.usage()

    # We only check 'start' option
    if args[0] == 'start':
        # Check if there are at least four arguments
        if len(args) < 5:
            print "Insufficient arguments for 'start' command."
            try:
                win32serviceutil.usage()
            except:
                pass
            start_usage()
            sys.exit(1)

        # Check if we can get to the working directory
        try:
            os.chdir(args[1])
        except:
            print "Cannot change to working directory."
            try:
                win32serviceutil.usage()
            except:
                pass
            start_usage()
            sys.exit(1)
           

class FSindexService(win32serviceutil.ServiceFramework):
    _svc_name_ = "FSindexService"
    _svc_display_name_ = "FSindex Python Service"
    
    def __init__(self, args):

        win32serviceutil.ServiceFramework.__init__(self, args)
        self.args = args

        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = file(os.path.join('c:\\', 'FSIndex_out.log'),'a')
        sys.stderr = file(os.path.join('c:\\', 'FSIndex_out.log'),'a')

    def SvcStop(self):

        self.ReportServiceStatus(win32service.SERVICE_STOP_PENDING)
        self.SrchS.terminate_flag = SearchServer.SIGNAL 
 
    def SvcDoRun(self):

        # Assume that arguments have been scanned and that workdir is a valid directory
        workdir, daemonid, port, indexfile = [str(unicode_str) for unicode_str in self.args[1:5]]

        # Now redirect to log files
        sys.stdout.close()
        sys.stderr.close()
        logfile = os.path.join(workdir, "FSsearchs_%s.log" % daemonid)
        sys.stdout = file(logfile, 'a')
        sys.stderr = file(logfile, 'a')

        self.ReportServiceStatus(win32service.SERVICE_START_PENDING)
        self.SrchS = SearchServer.SearchServer(daemonid, port, indexfile, control_slaves=False)
        self.ReportServiceStatus(win32service.SERVICE_RUNNING)
        self.SrchS.start()
        self.ReportServiceStatus(win32service.SERVICE_STOPPED)

if __name__=='__main__':
    check_start_command_line()
    win32serviceutil.HandleCommandLine(FSindexService)
