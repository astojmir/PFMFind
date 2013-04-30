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


import socket
from cStringIO import StringIO

from pfmfind.search.SearchServer import DIE, GET_INDEX_DATA, \
     GET_SERVERS, SEARCH, SCORE_DISTR, send_obj, receive_obj 

class SearchClient(object):
    def __init__(self):
        self.detach()

    def _master_data(self, opcode, address=None, data=0):
        if address == None:
            address = (self.host, self.port)
        try: 
            try:
                sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                sock.connect(address)
                send_obj(sock, (opcode,data))
                mdata = receive_obj(sock)
            finally:
                sock.close()
        except Exception, inst:
            return None
        return mdata
        
    def poll(self):
        if not self.attached: return
        if self.get_ix_data() == None:
            return False
        else:
            return True

    def get_ix_data(self, address=None):
        if not self.attached and address == None:
            return None
        return self._master_data(GET_INDEX_DATA, address)

    def get_servers(self, address=None):
        if not self.attached and address == None:
            return None
        return self._master_data(GET_SERVERS, address)

    def attach(self, host, port):
        if self.attached:
            raise RuntimeError,\
                  "Must be detached before attaching"

        ix_data = self.get_ix_data((host,port))
        if ix_data == None: return False

        self.host = host
        self.port = port
        self.ix_data = ix_data
        self.servers = self.get_servers((host,port))
        self.attached = True

        for ixd in ix_data:
            self.fragments += ixd['fragments']
            self.bins += ixd['bins']
            self.unique_fragments += ixd['unique_fragments']

        return True

    def detach(self):
        self.host = None
        self.port = None
        self.servers = None
        self.ix_data = None
        self.attached = False

        self.fragments = 0
        self.bins = 0
        self.unique_fragments = 0

    def kill(self):
        if not self.attached: return
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.connect((self.host, self.port))
        send_obj(sock, (0,0))
        sock.close()
        self.detach()

    def search(self, search_args):
        if not self.attached: return

        # Convert E-value to P-value 
        for j,srch in enumerate(search_args):
            srch.extend([None, self.fragments])

        # Search
        results = self._master_data(SEARCH, address=None, data=search_args)
        if results == None:
            return [None] * len(search_args)

        return results
    
    def index_data_str(self): 
        if self.ix_data == None:
            return ""

        fs = StringIO()

        fs.write("********* OVERALL SIZE *********\n")
        fs.write("Total fragments: %d\n" % self.fragments)
        fs.write("Total bins: %d\n\n" % self.bins)
        
        for j,ixdict in enumerate(self.ix_data):
            fs.write("********* INDEX #%d *********\n" % j)
            fs.write("***** Database Details *****\n")
            fs.write("Database name: %s\n" % ixdict['db_name'])
            fs.write("Full length: %d\n" % ixdict['db_length'])
            fs.write("Number of sequences: %d\n\n" % ixdict['db_no_seq'])
            fs.write("***** Index Details *****\n")
            fs.write("Index name: %s\n" % ixdict['index_name'])
            fs.write("Indexed Alphabet: %s\n" % ixdict['alphabet'])
            fs.write("Partitions:\n")
            ptable = ixdict['ptable']
            l = len(ptable)
            for i in range(0,l-2,2):
                fs.write((("%2.2d. " % i) + ptable[i] + "    "))
                fs.write((("%2.2d. " % (i+1)) + ptable[i+1] + "\n"))

            if l % 2 == 0:
                fs.write((("%2.2d. " % (l-2)) + ptable[l-2] + "    "))
            fs.write((("%2.2d. " % (l-1)) + ptable[l-1] + "\n"))

            fs.write("Number of bins: %d\n" % ixdict['bins'])
            fs.write("Number of indexed fragments: %d\n" % ixdict['fragments'])
            fs.write("Indexed fragment length: %d\n\n" % ixdict['indexed_fragment_length'])
            fs.write("\n\n")
        return fs.getvalue()
