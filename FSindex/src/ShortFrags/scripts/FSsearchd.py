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
Starts the FSsearch daemon that loads an instance of FSindex into memory
and waits for search requests.

SYNOPSIS:

    FSsearchd.py daemonid port workpath indexfile [pythonpath] 


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

"""


import sys, os, signal, os.path, string, resource, socket
from ShortFrags.Expt import SearchServer

__version__ = "$Revision: 1.1 $"
__date__ = "$Date$"
__author__ = "Aleksandar Stojmirovic"
__credits__ = "Chad J. Schroeder, ActiveState Programming Network"


# Based on Chad J. Schroeder's recipe: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/278731 

SrchS = None

def signal_handler(signum, frame):
    global SrchS
    SrchS.terminate_flag = SearchServer.SIGNAL

def create_daemon(daemonid, port, workpath, indexfile, pythonpath=None):
    """Disk And Execution MONitor (Daemon)

    Default daemon behaviors (they can be modified):
    1.) Ignore SIGHUP signals.
    2.) Default current working directory to the "/" directory.
    3.) Set the current file creation mode mask to 0.
    4.) Close all open files (0 to [SC_OPEN_MAX or 256]).
    5.) Redirect standard I/O streams to "/dev/null".

    Failed fork() calls will return a tuple: (errno, strerror).  This behavior
    can be modified to meet your program's needs.

    Resources:
    Advanced Programming in the Unix Environment: W. Richard Stevens
    Unix Network Programming (Volume 1): W. Richard Stevens
    http://www.erlenstar.demon.co.uk/unix/faq_2.html#SEC16
    """
    
    if pythonpath != None:
        sys.path.append(pythonpath)
    workpath = os.path.expanduser(workpath)

    try:
        # Fork a child process so the parent can exit.  This will return control
        # to the command line or shell.  This is required so that the new process
        # is guaranteed not to be a process group leader.  We have this guarantee
        # because the process GID of the parent is inherited by the child, but
        # the child gets a new PID, making it impossible for its PID to equal its
        # PGID.
        pid = os.fork()
    except OSError, e:
        return((e.errno, e.strerror))     # ERROR (return a tuple)

    if (pid == 0):       # The first child.

        # Next we call os.setsid() to become the session leader of this new
        # session.  The process also becomes the process group leader of the
        # new process group.  Since a controlling terminal is associated with a
        # session, and this new session has not yet acquired a controlling
        # terminal our process now has no controlling terminal.  This shouldn't
        # fail, since we're guaranteed that the child is not a process group
        # leader.
        os.setsid()

        # When the first child terminates, all processes in the second child
        # are sent a SIGHUP, so it's ignored.
        signal.signal(signal.SIGHUP, signal.SIG_IGN)

        try:
            # Fork a second child to prevent zombies.  Since the first child is
            # a session leader without a controlling terminal, it's possible for
            # it to acquire one by opening a terminal in the future.  This second
            # fork guarantees that the child is no longer a session leader, thus
            # preventing the daemon from ever acquiring a controlling terminal.
            pid = os.fork()        # Fork a second child.
        except OSError, e:
            return((e.errno, e.strerror))  # ERROR (return a tuple)

        if (pid == 0):      # The second child.
            # Ensure that the daemon doesn't keep any directory in use.  Failure
            # to do this could make a filesystem unmountable.
            os.chdir("/")
            # Give the child complete control over permissions.
            os.umask(0)
        else:
            pidfile = os.path.join(workpath, "FSsearchd_%s.pid" % daemonid)
            pid_fp = file(pidfile,'w')
            pid_fp.write("%d\n" % pid) # Write pid of the daemon
            pid_fp.write("%s\n" % socket.gethostname()) 
            pid_fp.write("%s\n" % string.join(sys.argv, ' ')) # Write the command line
            pid_fp.close()
            os._exit(0)     # Exit parent (the first child) of the second child.
    else:
        os._exit(0)         # Exit parent of the first child.



    # Change to working directory
    logfile = os.path.join(workpath, "FSsearchd_%s.log" % daemonid)
    os.chdir(workpath)

    # Close all open files.  Try the system configuration variable, SC_OPEN_MAX,
    # for the maximum number of open files to close.  If it doesn't exist, use
    # the default value (configurable).

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

    # Set resource limits
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
    SrchS = SearchServer.SearchServer(daemonid, port, indexfile)

    # Catch termination signal so we can write the log and kill all slaves
    signal.signal(signal.SIGTERM, signal_handler)
    
    SrchS.start()

    os._exit(0)
    
if __name__ == "__main__":
    if len(sys.argv) < 5:
        print __doc__
        sys.exit()
    args = tuple(sys.argv[1:])
    create_daemon(*args)
