import types, cPickle, threading, sys, os, re, socket, Queue
import time
from ShortFrags.Expt.index import FSIndex
from cStringIO import StringIO

from ShortFrags.Expt.matrix import SCORE, POSITIONAL
from ShortFrags.Expt.matrix import ScoreMatrix, ProfileMatrix
from ShortFrags.Expt.scoredistr import ScoreDistr


RNG_SRCH = 0
KNN_SRCH = 1
REL_SRCH = 2


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


class DispatchedSearch(threading.Thread):
    def __init__(self, host, port, query, queue):
        self._host = host
        self._port = port
        self._query = query
        self._queue = queue
        threading.Thread.__init__(self)

    def run(self):
        try:
            res1 = send_obj2(self._host, self._port, self._query)
        except Exception, inst:
            res1 = None
        self._queue.put(res1)


def now():
    return time.ctime(time.time())

class SearchServer(object):
    def __init__(self, servers):
        if len(servers) == 1:
            self.distr = False
        else:
            self.distr = True

        self.servers = servers
        self.sp_acc = re.compile('\((\w+)\)')

    def load_index(self):
        print "Loading Index at %s" % now()
        sys.stdout.flush()
        index_filename = self.servers[0][3]
        self.I = FSIndex(index_filename, [])

    def create_subservers(self):
        print "Creating subservers at %s" % now()
        sys.stdout.flush()
        self.queue = Queue.Queue(len(self.servers)-1)
        for i,srvr in enumerate(self.servers[1:]):
            fn = os.path.join(os.path.expanduser(srvr[2]), 'server%d.cfg' % i)
            fp = file(fn, 'w')
            fp.write('%s %d %s %s\n' % srvr)
            fp.close()
            command = 'rsh -n %s "FSsearchd.py %s >&/dev/null </dev/null &"'\
                      % (srvr[0], fn)
            os.system(command)

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

                if self.distr:
                    self.dispatch_request(clsock, opcode, data)
                else:
                    self.process_request(clsock, opcode, data)
            finally:
                clsock.close()
        except Exception, inst:
            print ("Failed request at %s: " % now()) , inst
            sys.stdout.flush()

    def dispatch_request(self, clsock, opcode, data):
        if opcode == DIE:
            print "Terminating all subservers at %s" % now() 
            sys.stdout.flush()
            for srvr in self.servers[1:]:
                try:
                    send_obj2(srvr[0], srvr[1], (opcode, data))
                except Exception, inst:
                    print "Response error: ", srvr, inst
                    sys.stdout.flush()
            print "Terminating at %s" % now()
            sys.stdout.flush()
            clsock.close()
            os._exit(0)
        elif opcode == GET_INDEX_DATA:
            results = []
            for srvr in self.servers[1:]:
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

            print "Started %d search requests at %s" % (len(data), now())
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

            chunk = (len(sd_data)//(len(self.servers)-1)) + 1
            a = 0
            b = chunk
            for j,srvr in enumerate(self.servers[1:]):
                DispatchedSearch(srvr[0], srvr[1],
                                 (SCORE_DISTR, ((a,b),sd_data)), self.queue).start()
                a += chunk
                b += chunk
                
            for j in xrange(len(self.servers)-1):
                rng, res = self.queue.get()
                for i in xrange(*rng):
                    data[i][6] = res[i]
                    

            # Use threads to produce the results simultaneously
            for srvr in self.servers[1:]:
                DispatchedSearch(srvr[0], srvr[1], (opcode, data), self.queue).start()

            # First search separately
            results = self.queue.get()

            # Must retrieve all results in order to clear the queue
            # Assume that all results are in order
            for j in xrange(len(self.servers)-2):
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

            # Tidy up if kNN search
            for i,HL in enumerate(results):
                if HL == None:
                    continue
                search_type = data[i][0]
                if search_type == KNN_SRCH:
                    k = data[i][4]  
                    HL.sort_by_distance()
                    d = HL[k-1].dist                    
                    while HL[k].dist == d:
                        k += 1
                    del(HL[k:])
            print "Processed %d search requests at %s" % (len(data), now())
            sys.stdout.flush()
        else:
            # Just do it
            self.process_request(clsock, opcode, data)

        send_obj(clsock, results)
        
    def process_request(self, clsock, opcode, data):
        if opcode == DIE:
            print "Terminating at %s" % now()
            sys.stdout.flush()
            clsock.close()
            os._exit(0)
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
            for srch in data:
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
                    HL=self.I.rng_srch(qseq, M, r0, conv_type)
                elif search_type == KNN_SRCH:
                    HL=self.I.kNN_srch(qseq, M, r0)
                elif search_type == REL_SRCH:
                    r1 = SD.cutoff(r0/totfrags)
                    HL=self.I.rng_srch(qseq, M, r1, conv_type)
                else:
                    results.append(None)
                    continue

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
            results = []
            rng = data[0]
            sd_data = data[1]
            for dargs in sd_data:
                if type(dargs) == types.TupleType and len(dargs) == 4:
                    results.append(ScoreDistr(*dargs))
                else:
                    results.append(None)
            results = (rng,results)
        send_obj(clsock, results)

            
    def server_main(self):
        print "Starting Main Loop at %s" % now()
        sys.stdout.flush()
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        host = ''
        port = self.servers[0][1]
        sock.bind((host, port))
        sock.listen(5)

        while 1:
            (clsock, address) = sock.accept()
            ct = threading.Thread(target=self.answer_request, args=(clsock,))
            ct.start()

    def start(self):
        print "Starting at %s" % now()
        sys.stdout.flush()
        if self.distr:
            self.create_subservers()
            self.server_main()
        else:
            self.load_index()
            self.server_main()



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
