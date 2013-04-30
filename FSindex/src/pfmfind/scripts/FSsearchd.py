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
Starts the FSsearch daemon that loads an instance of FSindex into memory
and waits for search requests.

SYNOPSIS:

    FSsearchd.py [--control-slaves] daemonid port workpath indexfile [pythonpath]
    FSsearchd.py -h|--help

DESCRIPTION:

This script starts a daemon that loads indexfile into memory and then
listens at a given port for search requests. The logfiles are written in
the workpath. The daemonid is an arbitrary identifier that is appended
to the logfiles written by each daemon. The optional argument pythonpath
can be used to add a directory to the default python path.

If the inxdexfile ends with .cfg it is treated as a configuration file
for slaves. Each line should be in the following format:

host port workpath indexfile [pythonpath [binpath]].

The fields should be separated by blanks.

OPTIONS:

--help (-h)             Shows this message.
--control-slaves (-c)   The server attempts to start its slaves if not already
                          started (using ssh). On shutdown, it sends the
                          termination signal to all slaves. Note that you need
                          to have passwordless ssh setup for this to work.
"""


import sys, os, signal, os.path, string, resource, socket
from pfmfind.Expt import SearchServer

__version__ = "$Revision: 1.1 $"
__date__ = "$Date$"
__author__ = "Aleksandar Stojmirovic"

# Based on the recipies from the Active State Programming Network Python Cookbook:
# http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/278731 (Chad J. Schroeder) and
# http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/66012 (Jurgen Hermann)

SrchS = None

def signal_handler(signum, frame):
    global SrchS
    SrchS.terminate_flag = SearchServer.SIGNAL

def create_daemon(daemonid, port, workpath, indexfile, pythonpath=None, control_slaves = False):
    """
    Creates a daemon to run as PFMFind FSIndex server.
    """

    if pythonpath != None:
        sys.path.append(pythonpath)
    workpath = os.path.expanduser(workpath)
    logfile = os.path.join(workpath, "FSsearchd_%s.log" % daemonid)
    pidfile = os.path.join(workpath, "FSsearchd_%s.pid" % daemonid)

    # First check if a daemon is already running, exit if so

    try:
        # Get pid from the pid file
        fp = file(pidfile, 'r')
        pid = int(fp.next().strip())
        fp.close()
        os.kill(pid, 0)
        sys.stderr.write("Another instance of FSsearchd running - kill it first.\n")
        os._exit(1)
    except:
        # Could not open file, could not find pid or the process itself - this is good
        pass

    try: # First fork
        pid = os.fork()
    except OSError, e:
        raise Exception, "%d (%s)" % (e.strerror, e.errno)

    if (pid == 0):
        # We are now in the first child.
        os.setsid()
        os.chdir("/")
        os.umask(0)
        signal.signal(signal.SIGHUP, signal.SIG_IGN)

        try:  # Fork a first child.
            pid = os.fork()
        except OSError, e:
            raise Exception, "%d (%s)" % (e.strerror, e.errno)

        if (pid > 0):
            # The parent (the first child) writes the pid of the daemon,
            # the host and the command line, then exits
            pid_fp = file(pidfile,'w')
            pid_fp.write("%d\n" % pid)
            pid_fp.write("%s\n" % socket.gethostname())
            pid_fp.write("%s\n" % string.join(sys.argv, ' '))
            pid_fp.close()
            os._exit(0)     # Exit parent of the second child.
    else:
        os._exit(0)         # Exit parent of the first child.


    # We are now in the second child.

    # Change to working directory
    os.chdir(workpath)

    # Close all open files.  Try the system configuration variable, SC_OPEN_MAX,
    # for the maximum number of open files to close.  If it doesn't exist, use
    # the default value (256).

    try:
        maxfd = os.sysconf("SC_OPEN_MAX")
    except (AttributeError, ValueError):
        maxfd = 256       # default maximum

    for fd in range(0, maxfd):
        try:
            os.close(fd)
        except OSError:   # ERROR (ignore)
            pass

    # Redirect the standard file descriptors to /dev/null and logfiles.
    os.open("/dev/null", os.O_RDONLY)                  # standard input  (0)
    os.open(logfile, os.O_CREAT|os.O_APPEND|os.O_RDWR) # standard output (1)
    os.open(logfile, os.O_CREAT|os.O_APPEND|os.O_RDWR) # standard error  (2)

    # Set resource limits - FSIndex takes a lot of memory
    if 'RLIMIT_DATA' in dir(resource):
        lim = resource.getrlimit(resource.RLIMIT_DATA)
        resource.setrlimit(resource.RLIMIT_DATA, (lim[1],lim[1]))
    if 'RLIMIT_RSS' in dir(resource):
        lim = resource.getrlimit(resource.RLIMIT_RSS)
        resource.setrlimit(resource.RLIMIT_RSS, (lim[1],lim[1]))
    if 'RLIMIT_MEMLOCK' in dir(resource):
        lim = resource.getrlimit(resource.RLIMIT_MEMLOCK)
        resource.setrlimit(resource.RLIMIT_MEMLOCK, (lim[1],lim[1]))

    # Start the SearchServer
    global SrchS
    SrchS = SearchServer.SearchServer(daemonid, port, indexfile, control_slaves)

    # Catch termination signal so we can write the log and kill all slaves
    signal.signal(signal.SIGTERM, signal_handler)

    SrchS.start()

    # Exit when the server is stopped
    os._exit(0)

if __name__ == "__main__":

    import getopt

    control_slaves = False
    pythonpath = None

    options = 'hc'
    long_options = ['help', 'control-slaves']

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
        elif o in ('-c', '--control-slaves'):
            control_slaves = True

    if len(args) < 4:
        print __doc__
        sys.exit(2)

    daemonid, port, workpath, indexfile = args[:4]

    if len(args) >= 5:
        pythonpath = args[4]

    create_daemon(daemonid, port, workpath, indexfile, pythonpath, control_slaves)
