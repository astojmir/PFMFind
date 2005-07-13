#! /usr/local/bin/python

import socket, sys, time, string
from ShortFrags.Expt.SearchServer import *

def now():
    return time.ctime(time.time())

def kill_all_slaves(serverfile):
    print "Terminating all subservers at %s" % now() 
    sys.stdout.flush()

    # Read the configuration file
    # File format: host port workpath indexfile [pythonpath [binpath]]
    fp = file(serverfile, 'r')
    servers = []
    for line in fp:
        sp_line = string.split(line)
        sp_line[1] = int(sp_line[1])
        servers.append(sp_line)
    fp.close()

    for srvr in servers:
        try:
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.connect((srvr[0], srvr[1]))
            send_obj(sock, (0, 0))
        except Exception, inst:
            print "Response error: ", srvr, inst
            sys.stdout.flush()
        print "Terminating at %s" % now()
        sys.stdout.flush()

if __name__ == "__main__":
    kill_all_slaves(sys.argv[1])
