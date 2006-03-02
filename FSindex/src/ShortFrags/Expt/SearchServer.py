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


import types, cPickle, sys, os, re, socket, Queue, select
import time, os.path, string
from ShortFrags.Expt.index import FSIndex
from cStringIO import StringIO
from threading import Thread
from errno import EINTR

from ShortFrags.Expt.matrix import SCORE, POSITIONAL
from ShortFrags.Expt.matrix import ScoreMatrix, ProfileMatrix
from ShortFrags.Expt.scoredistr import ScoreDistr


RNG_SRCH = 0
KNN_SRCH = 1
REL_SRCH = 2

MAX_HITS = 1500


# ***********************************************************************
# CLIENT-SERVER COMMUNICATIONS PROTOCOL
#  
# Client sends requests to server as pickled Python objects, server
# answers back the same way. A single request is allowed through
# a communications channel.
#
# A request is a 2-tuple of the form (opcode, data) where data depends
# on opcode. Server returns what it needs to, or None object if it
# cannot satisfy the request.
#
# OPCODES:
#
# DIE -            the server and any of its subservers exit, returns
#                  nothing. 
# GET_INDEX_DATA - returns a list (for each of subservers) of index
#                  descriptions. The number of descriptions depends on
#                  topology (i.e. a subserver can have more subservers
#                  ...). None is returned if some of subservers cannot
#                  be reached.
# GET_SERVERS -    returns the servers list for the current server.
# SEARCH -         Data is in the form of a list of 7-tuples, each tuple
#                  denoting a search. The tuple entries are:
#                  0 - search_type: RNG_SRCH or KNN_SRCH
#                  1 - qseq: Query sequence
#                  2 - matrix: Biopython object convertible into a matrix
#                  3 - matrix_type: SCORE or  POSITIONAL
#                  4 - r0: cutoff value, depends on search_type
#                  5 - conv_type: matrix conversion type
#                  6 - SD: score distribution. If None, will calculate it.
#                  7 - totfrags: Total number of fragments in database
#
#                  For each search, a tuple (id, HL) where HL is hit
#                  list or None (if search cannot be done) is returned.
#                  The tuples are packed into a list.
# SCORE_DISTR -    Calculate score distributions for a list of queries.
# ***********************************************************************


# OPCODES
DIE = 0
GET_INDEX_DATA = 1
GET_SERVERS = 2
SEARCH = 3
SCORE_DISTR = 4

def send_obj(sock, obj):
    msg = cPickle.dumps(obj, 2)   
    totalsent = 0
    while totalsent < len(msg):
        sent = sock.send(msg[totalsent:])
        if sent == 0:
            raise RuntimeError, "Socket connection broken"
        totalsent += sent
    sock.shutdown(1)

def send_obj2(host, port, obj):
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        sock.connect((host, port))
        send_obj(sock, obj)
        result = receive_obj(sock)
    finally:
        sock.close()
    return result

def receive_obj(sock):
    msg_str = StringIO()
    while 1:
        chunk = sock.recv(4096)
        if chunk == '':
            break
        msg_str.write(chunk)
    return cPickle.loads(msg_str.getvalue())


def parse_slaves_config(serverfile):
    # File format: host port workpath indexfile [pythonpath [binpath]]
    fp = file(serverfile, 'r')
    servers = []
    for line in fp:
        sp_line = string.split(line)
        sp_line[1] = int(sp_line[1])
        servers.append(sp_line)
    fp.close()
    return servers


def terminate_slaves(servers):
    if not servers:
        return
    print "Terminating all subservers at %s" % now()
    sys.stdout.flush()
    for srvr in servers:
        try:
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.connect((srvr[0], srvr[1]))
            send_obj(sock, (0, 0))
        except Exception, inst:
            print "Response error: ", srvr, inst
            sys.stdout.flush()





class DispatchedSearch(Thread):
    def __init__(self, host, port, query, queue):
        self._host = host
        self._port = port
        self._query = query
        self._queue = queue
        Thread.__init__(self)

    def run(self):
        try:
            res1 = send_obj2(self._host, self._port, self._query)
        except Exception, inst:
            res1 = None
        self._queue.put(res1)


def now():
    return time.ctime(time.time())


# Termination flag values:
RUN = 0           # keep running
REMOTE = 1        # remote termination - DIE request
SIGNAL = 2        # received termination signal
NO_SLAVES = 3     # could not reach a slave
NO_INDEX = 4      # could not load index
OTHER_ERROR = 99  # anything not numbered above

class SearchServer(object):
    def __init__(self, daemonid, port, indexfile, control_slaves=True):

        self.daemonid = daemonid
        self.port = int(port)
        self.indexfile = indexfile

        self.ct = None

        self.terminate_flag = RUN
        self.control_slaves = control_slaves

    def start(self):
        """
        Starts the daemon main loop.
        """

        print "Starting at %s" % now()
        sys.stdout.flush()

        # Check if indexfile is in fact a config file
        if os.path.splitext(self.indexfile)[1] == '.cfg':
            self.create_subservers(self.indexfile)
        else:
            self.servers = []
            self.load_index(self.indexfile)
        
        if not self.terminate_flag:
            self.sp_acc = re.compile('\((\w+)\)')

            print "Starting Main Loop at %s" % now()
            sys.stdout.flush()
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.bind(('', self.port))
            sock.listen(5)

        while not self.terminate_flag:
            try:
                r, w, e = select.select([sock],[],[], 2.0)
            except select.error, inst:
                if inst[0] != EINTR:
                    print "Select error in the main loop:" , inst.__class__, inst.__str__()
                    sys.stdout.flush()
                    self.terminate_flag = OTHER_ERROR
                r = []

            if r:
                (clsock, address) = sock.accept()
                self.ct = Thread(target=self.answer_request, args=(clsock,))
                self.ct.start()
                # Should not run any requests concurrently - only one additional thread
                # TO DO: make sure this thread always finishes.
                self.ct.join()
                self.ct = None

        self.terminate(terminate_subservers=self.control_slaves)

    def terminate(self, terminate_subservers=True):
        """
        Terminates cleanly (writes to log, terminates subservers).
        """

        if self.terminate_flag == REMOTE:
            print "Exiting due to DIE request."
        elif self.terminate_flag == SIGNAL:
            print "Exiting due to termination signal."
        elif self.terminate_flag == NO_SLAVES:
            print "Exiting due to being unable to reach all slaves."
        elif self.terminate_flag == NO_INDEX:
            print "Exiting due to being unable to load index."
            

            
        sys.stdout.flush()

        # Wait for working thread to finish
        if self.ct:
            print "Waiting for working thread at %s" % now()
            self.ct.join()

        if terminate_subservers:
            terminate_slaves(self.servers)
        print "Terminating at %s" % now()
        sys.stdout.flush()


    def load_index(self, indexfile):
        print "Loading Index at %s" % now()
        sys.stdout.flush()
        try:
            self.I = FSIndex(indexfile, [])
        except Exception, inst:
            print "Error loading index:" , inst.__class__, inst.__str__()
            sys.stdout.flush()
            self.terminate_flag = NO_INDEX

    def create_subservers(self, serverfile):
        print "Creating subservers at %s" % now()
        sys.stdout.flush()

        # Read the configuration file
        self.servers = parse_slaves_config(serverfile)

        self.queue = Queue.Queue(len(self.servers))

        for i,srvr in enumerate(self.servers):
            host, port, workpath, indexfile = srvr[:4]
            pythonpath = binpath = ''            
            if len(srvr) > 4:
                pythonpath = srvr[4]
                if len(srvr) > 5:
                    binpath = srvr[5]
            try: # Try to see if a daemon is already running
                send_obj2(host, port, (1,0))
                print "Found FSsearch slave %2.2d running on %s:%d" % (i, host, port)
                sys.stdout.flush()
            except socket.error, inst:
                if self.control_slaves: # Start the daemon 
                    print "Starting FSsearch slave %2.2d on %s:%d using ssh." % (i, host, port)
                    sys.stdout.flush()
                    daemon_full = os.path.join(binpath, "FSsearchd.py")
                    daemon_args = '%s_s%2.2d %d %s %s %s' % (self.daemonid, i, port, workpath,
                                                             indexfile, pythonpath)
                    command = 'ssh -x %s "%s %s >&/dev/null </dev/null &"'\
                              % (host, daemon_full, daemon_args)
                    os.system(command)

                    # (TO DO:) We should probably loop here (or elsewhere) until all servers are
                    # on line. On the other hand, clients check that before each search.

                else: # We just exit
                    print "Could not find FSsearch slave %2.2d running on %s:%d" % (i, host, port)
                    sys.stdout.flush()
                    self.terminate_flag = NO_SLAVES

    def check_request(self, clsock):
        obj = receive_obj(clsock)
        if type(obj) != types.TupleType and len(obj) != 2:
            raise RuntimeError, "Wrong type recieved"

        opcode = obj[0]
        data = obj[1]
                
        if opcode > 4:
            raise RuntimeError, "Invalid opcode recieved"

        return opcode, data

    def answer_request(self, clsock):
        try:
            try:
                opcode, data = self.check_request(clsock)
                print "Answering request %d at %s: " % (opcode, now())
                sys.stdout.flush()

                if len(self.servers):
                    self.dispatch_request(clsock, opcode, data)
                else:
                    self.process_request(clsock, opcode, data)
            finally:
                clsock.close()
        except Exception, inst:
            print ("Failed request at %s: " % now()) , inst.__class__, inst.__str__()
            sys.stdout.flush()

    def dispatch_request(self, clsock, opcode, data):
        if opcode == DIE:
            self.terminate_flag = REMOTE
            return
        elif opcode == GET_INDEX_DATA:
            results = []
            for srvr in self.servers:
                try:
                    res1 = send_obj2(srvr[0], srvr[1], (opcode, data))
                    results.extend(res1)
                except Exception, inst:
                    print "Response error: ", srvr, inst
                    sys.stdout.flush()
                    raise RuntimeError, "Subserver Failure"
        elif opcode == GET_SERVERS:
            results = self.servers
        elif opcode == SEARCH:
            if type(data) != types.ListType:
                raise RuntimeError,\
                      "Wrong search request type recieved"

            print "Started %d distribution requests at %s" % (len(data), now())
            sys.stdout.flush()

            # Dispatch calculation of SD
            sd_data = []
            for srch in data:
                if type(srch) != types.ListType\
                       or len(srch) != 8:
                    sd_data.append(None)
                else:
                    qseq = srch[1]
                    matrix = srch[2]
                    matrix_type = srch[3]
                    conv_type = srch[5]
                    sd_data.append((matrix, matrix_type, conv_type, qseq))

            chunk = len(sd_data)//len(self.servers)
            if chunk > 0:
                for j,srvr in enumerate(self.servers):
                    a = j*chunk
                    b = a+chunk
                    DispatchedSearch(srvr[0], srvr[1],
                                     (SCORE_DISTR, ((a,b), sd_data[a:b])),
                                     self.queue).start()
            else:
                j = 0

            a = len(self.servers)*chunk
            rng, res = self._process_distributions(((a,len(sd_data)),sd_data[a:]))
            for j,i in enumerate(xrange(*rng)):
                data[i][6] = res[j]
                
            if chunk > 0:
                for srvr in self.servers:
                    rng, res = self.queue.get()
                    for j,i in enumerate(xrange(*rng)):
                        data[i][6] = res[j]
                    

            print "Started %d search requests at %s" % (len(data), now())
            sys.stdout.flush()
            # Use threads to produce the results simultaneously
            for srvr in self.servers:
                DispatchedSearch(srvr[0], srvr[1], (opcode, data), self.queue).start()

            # First search separately
            results = self.queue.get()

            # Must retrieve all results in order to clear the queue
            # Assume that all results are in order
            for j in xrange(len(self.servers)-1):
                res1 = self.queue.get()
                if res1 == None: results = None
                if results == None: continue
                for i,HL in enumerate(res1):
                    if HL == None or results[i] == None:
                        results[i] = None
                        continue
                    results[i].extend(HL)
                    results[i].frags_visited += HL.frags_visited
                    results[i].frags_hit += HL.frags_hit
                    results[i].bins_visited += HL.bins_visited
                    results[i].bins_hit += HL.bins_hit
                    results[i].unique_frags_visited += HL.unique_frags_visited
                    results[i].unique_frags_hit += HL.unique_frags_hit
                
            # Check that the queue is empty
            if not self.queue.empty():
                print "Queue not empty"
                sys.stdout.flush()

            if results == None:
                raise RuntimeError, "Subserver Failure"

            print "Collecting kNN searches at %s" % now()
            sys.stdout.flush()

            # Tidy up if kNN search
            for i,HL in enumerate(results):
                if HL == None:
                    continue
                search_type = data[i][0]
                if search_type == KNN_SRCH:
                    k = data[i][4]  
                    HL.sort_by_distance()
                    d = HL[k-1].dist                    
                    if HL[-1].dist > d:
                        while HL[k].dist == d:
                            k += 1
                        del(HL[k:])
            print "Processed %d search requests at %s" % (len(data), now())
            sys.stdout.flush()
        else:
            results = self._process_distributions(data)

        send_obj(clsock, results)

    def _process_distributions(self, data):
        results = []
        rng = data[0]
        sd_data = data[1]
        for dargs in sd_data:
            if type(dargs) == types.TupleType and len(dargs) == 4:
                results.append(ScoreDistr(*dargs))
            else:
                results.append(None)
        results = (rng,results)
        print "Processed %d distributions at %s" % (len(sd_data), now())
        sys.stdout.flush() 
        return results
        
    def process_request(self, clsock, opcode, data):
        if opcode == DIE:
            self.terminate_flag = REMOTE
            return
        elif opcode == GET_INDEX_DATA:
            results = [self.I.ix_data]
        elif opcode == GET_SERVERS:
            results = self.servers
        elif opcode == SEARCH:
            if type(data) != types.ListType:
                raise RuntimeError,\
                      "Wrong search request type recieved"

            print "Started %d search requests at %s" % (len(data), now())
            sys.stdout.flush()

            results = []
            for i, srch in enumerate(data):
                if type(srch) != types.ListType\
                       or len(srch) != 8:
                    results.append(None)
                    continue
                
                search_type = srch[0]
                qseq = srch[1]
                matrix = srch[2]
                matrix_type = srch[3]
                r0 = srch[4]
                conv_type = srch[5]
                SD = srch[6]
                totfrags = srch[7]
                
                if matrix_type == SCORE:
                    M = ScoreMatrix(matrix)
                elif  matrix_type == POSITIONAL:
                    M = ProfileMatrix(matrix.pssm)
                else:
                    results.append(None)
                    continue
                
                if SD == None:
                    SD = srch[6] = ScoreDistr(matrix, matrix_type, conv_type, qseq)
                
                if conv_type:
                    M.conv_type = conv_type
                    M = M.matrix_conv()

                if search_type == RNG_SRCH:
                    HL=self.I.rng_srch(qseq, M, r0)
                elif search_type == KNN_SRCH:
                    HL=self.I.kNN_srch(qseq, M, r0)
                elif search_type == REL_SRCH:
                    r1 = SD.cutoff(r0/totfrags)
                    HL=self.I.rng_srch(qseq, M, r1)
                else:
                    results.append(None)
                    continue

                # Truncate too large sets of results.
                if len(HL) > MAX_HITS:
                    HL.sort_by_distance_only()
                    del(HL[MAX_HITS:])
                    print "Truncated hit list %d" % i
                    sys.stdout.flush()


                # Now add sequence and accession to the hit
                # At some stage this should be implemented in
                # C Also note that this only works for
                # SwissProt/TrEMBL 

                for hit in HL:
                    # Set accession
                    hit.sequence = self.I.db.get_frag(hit.seq_id, hit.seq_from, hit.seq_to)                  
                    defline = self.I.db.get_def(hit.seq_id)
                    m = self.sp_acc.search(defline) 
                    hit.accession = m.groups()[0]
                    del(hit.seq_id)

                    # Set E-value and p-value
                    if conv_type == 0: 
                        score = hit.sim
                    else:
                        score = hit.dist
                    hit.pvalue = SD.pvalue(score)
                    hit.Evalue = totfrags * hit.pvalue
                   

                results.append(HL)
            print "Processed %d search requests at %s" % (len(data), now())
            sys.stdout.flush()
        else:
            results = self._process_distributions(data)
        send_obj(clsock, results)
