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
Kills all subservers.

"""

## import socket, sys, time, string

from pfmfind.Expt.SearchServer import parse_slaves_config
from pfmfind.Expt.SearchServer import terminate_slaves


## def now():
##     return time.ctime(time.time())

## def kill_all_slaves(serverfile):
##     print "Terminating all subservers at %s" % now() 
##     sys.stdout.flush()

##     # Read the configuration file
##     # File format: host port workpath indexfile [pythonpath [binpath]]
##     fp = file(serverfile, 'r')
##     servers = []
##     for line in fp:
##         sp_line = string.split(line)
##         sp_line[1] = int(sp_line[1])
##         servers.append(sp_line)
##     fp.close()

##     for srvr in servers:
##         try:
##             sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
##             sock.connect((srvr[0], srvr[1]))
##             send_obj(sock, (0, 0))
##         except Exception, inst:
##             print "Response error: ", srvr, inst
##             sys.stdout.flush()
##         print "Terminating at %s" % now()
##         sys.stdout.flush()

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print __doc__
        sys.exit()

    terminate_slaves(parse_slaves_config(sys.argv[1]))    
