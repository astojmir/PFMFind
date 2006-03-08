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



import socket, sys
from ShortFrags.Expt.SearchClient import SearchClient, send_obj
from ShortFrags.Expt.matrix import SCORE, POSITIONAL
from Bio.SubsMat import MatrixInfo


HOST = 'localhost'
PORT = 50007

matrix = MatrixInfo.blosum62

obj1 = (0,
        'YPQPQPI',
        matrix,
        SCORE,
        22,
        0,
        )

obj2 = (1,
        'IVLFFIV',
        matrix,
        SCORE,
        214,
        0,
        )




if __name__ == "__main__":
    test_type = int(sys.argv[1])

    if test_type == 0:
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.connect((HOST, PORT))
        send_obj(sock, (0,0))
        sock.close()
        sys.exit(0)

    SC = SearchClient()
    flag = SC.attach(HOST, PORT)

    if not flag:
        print "Could not attach!"
        sys.exit(1)

    if test_type == 1:
        # GET_INDEX_DATA test
        print SC.get_ix_data()
    elif test_type == 2: 
        # GET_SERVERS test
        print SC.get_servers()
    else:
        # SEARCH test
        result = SC.search([obj1,obj2])
        print len(result)
        #print result[0]
        print len(result[1])
        #print result[0][1][0].__dict__
        #print result[0][1].__dict__

        del(result[0][10:])
        print result[0]
        
