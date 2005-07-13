#!/usr/bin/env python


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
        
