"""
Provides the code to manipulate a full fragment experiment.

B{Classes:}
    - FullExpt
        - Provides the functions to load an index and 
        FASTA database, to run a search, display and save the results 
        of that experiment. 

B{Functions:}
    - print_matrix(matrix, conv_type, scale=1.0, weight_type=None, 
    dirichlet_type=None, FE=None, coords=None)
        - Return the distance matrix to be printed.
    - print_align(weight_type, dirichlet_type, FE, coords)
        - Return a printable string with the alignment information for 
        an instance of FullExpt.
    - load_FE(filename)
        - Load a saved experiment file.
    
B{Exceptions:}
    - IOError
"""
import os
import os.path
import types
import sys
import md5
from ShortFrags import FS
from hit_list import *
from index import IndexedDb
from DirichletMix import DirichletMix, freq_counts, henikoff_weights
from DirichletMix import BKGRND_PROBS as bg_dict
from DirichletInfo import get_mix
import threading
import cPickle
import gzip
from cStringIO import StringIO

from Bio import SubsMat
from Bio.SubsMat import MatrixInfo
from Bio.SubsMat.FreqTable import FreqTable, FREQ
from Bio.Alphabet import IUPAC
from Bio.Align import Generic, AlignInfo


RNG_SRCH = 0
KNN_SRCH = 1
REL_SRCH = 2

MATRIX_CTYPE = {'None': 0, 'Quasi': FS.QUASI,
                'Avg': FS.AVG, 'Max': FS.MAX}
		
matrix_cache = {}

def _get_converted_matrix(matrix, conv_type, scale=1.0, weight_type=None,
                         dirichlet_type=None, FE=None, coords=None):
    PM = None

    if weight_type != None: # Iterated search
        l = coords[0]
        f = coords[1]
        i = coords[2]

        # Assuming that the iteration is correct
        HL = FE.get_iters(l,f)[i-1]
        seqs = HL.get_seqs()
        
        if conv_type == -1:
            conv_type = 'None'

        # Calculate sequence weights
        bcounts = freq_counts(seqs, "ACDEFGHIKLMNPQRSTVWY")
        if weight_type == 'None':
            weights = [1.0]*len(seqs)
        elif weight_type == 'Henikoff':
            weights = henikoff_weights(seqs, "ACDEFGHIKLMNPQRSTVWY", bcounts)

        if matrix == 'PSSM': 
            DM = get_mix(dirichlet_type)
            if weight_type == 'Henikoff':
                bcounts = DM.block_counts(seqs, weights)
            bprobs = DM.block_probs(bcounts)
            bkgrnd = DM.aa_vector(bg_dict)
            PM = DM.block2pssm(DM.block_log_odds(bprobs, bkgrnd, scale), seqs[0])
            M0 = FS.Smatrix_pssm(PM.pssm)

        else: # Score matrix 
            align = Generic.Alignment(IUPAC.protein)
            for i in range(len(seqs)):
                align.add_sequence("Seq #%d" % i, seqs[i], weight=weights[i])
            summary_align = AlignInfo.SummaryInfo(align)

            # Must get expected frequencies from our own (i.e. whole
            # database) background frequencies. Otherwise they would be
            # derived from the alignment which wouldn't be good if we
            # have a small sample
            ftab = FreqTable(bg_dict, FREQ)
            arm = SubsMat.SeqMat(summary_align.replacement_dictionary())
            lom = SubsMat.make_log_odds_matrix(arm, ftab, factor=scale,
                                               round_digit=0, keep_nd=0)
            M0 = FS.Smatrix(lom)

    elif matrix in matrix_cache:
        M0 = matrix_cache[matrix]
    else:
        if matrix in MatrixInfo.__dict__:
            M0 = FS.Smatrix(MatrixInfo.__dict__[matrix])
            matrix_cache[matrix] = M0
        else:
            return ""
    ctype = MATRIX_CTYPE[conv_type]
    if ctype:
        M0.conv_type = ctype
        M1 = M0.matrix_conv('')
    else:
        M1 = M0
    return M1, ctype, PM

def print_matrix(matrix, conv_type, scale=1.0, weight_type=None,
                 dirichlet_type=None, FE=None, coords=None):
    """
    Return the score matrix to be printed.
    
    The last four arguments need only be specified if matrix == 'PSSM'. 
    
    @param matrix: Matrix name. Could be any of the matrices from the 
    Biopython MatrixInfo or 'PSSM'.
    @param conv_type: Matrix conversion type. Possible values are the 
    keys of the dictionary B{MATRIX_CTYPE}.
    @param scale: Scaling factor for matrix values.
    @param weight_type: Type of sequence weights (None or 'Henikoff').
    @param dirichlet_type: Name of Dirichlet mixture.  
    Must be contained in B{NAMES} list in DirichletMix.py.
    @param FE: A loaded experiment.
    @param coords: The triple (l,f,i) describing current coordinates in 
    the experiment. The fragment hits associated with (l,f,i-1) are used 
    to compute the iterative matrices. 
    
    """		 
    M, ctype, PM = _get_converted_matrix(matrix, conv_type, scale, weight_type,
                                        dirichlet_type, FE, coords)
    return M.__str__()

def print_align(weight_type, dirichlet_type, FE, coords):
    """
    Return a printable string with the multiple alignment information for 
    the iteration previous to the one given by coords = (l,f,i).
    
    @param weight_type: Type of sequence weights (None or 'Henikoff').
    @param dirichlet_type: Name of Dirichlet mixture.  
    Must be contained in B{NAMES} list in DirichletMix.py.
    @param FE: A loaded experiment.
    @param coords: The triple (l,f,i) describing current coordinates in the 
    experiment. The fragment hits associated with (l,f,i-1) are used to 
    compute the iterative matrices. 
    
    @return:  A printable string providing the alignment information.
    """
    if weight_type == None: return ""

    l = coords[0]
    f = coords[1]
    i = coords[2]

    # Assuming that the iteration is correct
    HL = FE.get_iters(l,f)[i-1]
    seqs = HL.get_seqs()
    deflines = HL.get_deflines()
    DM = get_mix(dirichlet_type)

    bcounts = DM.block_counts(seqs)
    if weight_type == 'None':
        weights = [1.0]*len(seqs)
        wcounts = bcounts
    elif weight_type == 'Henikoff':
        weights = henikoff_weights(seqs, "ACDEFGHIKLMNPQRSTVWY", bcounts)
        wcounts = DM.block_counts(seqs, weights)

    file_str = StringIO()
    file_str.write('***** ALIGNMENT *****\n')
    for i in range(len(seqs)):
        file_str.write('%8.4f %s %s\n' % (weights[i], seqs[i], deflines[i]))
    file_str.write('\n***** COUNTS *****\n')
    file_str.write(DM.print_block_data(bcounts))
    file_str.write('\n***** WEIGHTED COUNTS *****\n')
    file_str.write(DM.print_block_data(wcounts, 5, 1, 'float'))
    file_str.write('\n***** DIRICHLET MIXTURE PROBABILITIES *****\n')
    bprobs = DM.block_probs(wcounts)
    file_str.write(DM.print_block_data(bprobs, 6, 4, 'float'))
    #file_str.write('\n***** INFORMATION CONTENT *****\n')
    #file_str.write(DM.block2pssm(wcounts, seqs[0]).__str__())
    return file_str.getvalue()

def _filter_ext(s, ext):
    """
    Check to make sure that the file extension is 
    correct.  
    
    The file path and extension are passed
    into the function.  The ext passed into the 
    function is compared against the one in the file
    path and the result of that comparison is returned.
    
    @param s: The file path.
    @param ext: The file extension.
    """
    from os.path import splitext
    fsplit = splitext(s)
    return fsplit[1] == ext

def load_FE(filename):
    """
    Load a saved experiment.
    
    @param filename: The name of the experiment file to load.
    
    @return:  Return the opened experiment file.
    """
    try:
        fp = gzip.GzipFile(filename, 'rb')
        buffer = ""
        while 1:
            data = fp.read()
            if data == "":
                break
            buffer += data
            FE = cPickle.loads(buffer)
    except IOError:
        fp.close()
        raise IOError, 'Could not load the file %s.' % filename
    
    fp.close()
    return FE

class FullExpt:
    """
    Provides the functions to load an index and 
    FASTA database, to run a search, and save the results 
    of that experiment. 
    
    An instance should either by created by using load_FE or by calling 
    the constructor AND set_params. In any case  load_fastadb or load_index 
    need to be called for it to be useful.
    
    B{Exceptions:}
        - IOError
    """

    def __init__(self):
        """
	Constructor, initializes data.  To load an 
	already saved instance of FE load_FE must be 
	called. 
	"""
        self.mcache = {}
        self.Idata = None
        self.index_filename = None
        self.fasta_filename = None
        self.desc_filename = None
        self.fasta_md5sum = None
        self.query_seq = None
        
    def set_params(self, qseq, qdef='', qrange=None, lrange=None):
        """
	Set the parameters for the full experiment and initialize 
	the experiment cube.
	
	@param qseq:  The query sequence.
	@param qdef:  The description of the query sequence.
	@param qrange:  Range of the query.  If None the default
	value is used.  The default value is the same as the 
	length of the query sequence.
	@param lrange:  Range of the length.  If None the 
	default value of 9 is used.  
	"""
        self.query_seq = qseq
        self.query_def = qdef
        if qrange == None:
            qrange = (0,len(qseq))
        self.query_range = qrange
        self.rqseq = qseq[qrange[0]:qrange[1]]
        if lrange == None:
            lrange = (6,15)
        self.length_range = lrange

        # Initialise the experiment cube
        lsize = lrange[1] - lrange[0]
        fsize = qrange[1] - qrange[0]
        self.E = [[[] for i in range(fsize+1-j-lrange[0])]
                  for j in range(lsize)]

    def load_fastadb(self, filename):
        """
        Load the FASTA database.  The function initializes 
        B{Idata} by instantiating an instance of B{IndexedDb}.
        The database is then loaded by calling the B{load_fasta_db} 
        function in B{IndexedDb}.
       
        @param filename:  The name of the FASTA database to be loaded. 
        """
        if self.Idata == None:
            self.Idata = IndexedDb()
        self.index_filename = None
        self.fasta_filename = None
        self.fasta_md5sum = self.Idata.load_fasta_db(filename,
                                                     self.fasta_md5sum)
        self.fasta_filename = filename
        
    def load_index(self, filename):
        """
        Load the index.  The function initializes B{Idata} by 
        instantiating an instance of B{IndexedDb}.  The index is 
        then loaded by calling the B{load_index} function in 
        B{IndexedDb}.
       
        @param filename: The name of the index to be loaded.
        """
        if self.Idata == None:
            self.Idata = IndexedDb()
        self.index_filename = None
        self.fasta_filename = None
        self.fasta_md5sum = self.Idata.load_index(filename,
                                                  self.fasta_md5sum)
        self.index_filename = filename
       
    def load_descriptions(self, filename):
        """
        Load the descriptions.  The function initializes B{Idata} by 
        instantiating an instance of B{IndexedDB}. The descriptions
        are then loaded by calling the B{load_descriptions} function
        in B{IndexedDb}.
       
        @param filename: The name of the file containing the descriptions
        to be loaded.
        """
        if self.Idata == None:
            self.Idata = IndexedDb()
        self.desc_filename = None
        self.Idata.load_descriptions(filename)
        self.desc_filename = filename

    def save(self, filename):
        """
        Save the experiment under the given filename.
        Raise an exception if the file cannot
        be saved.
       
        @param filename: The name the file will be saved under.
        """
        tmp = self.Idata
        self.Idata = None
        try:
            fp = gzip.GzipFile(filename, 'wb')
            fp.write(cPickle.dumps(self, 2))
        except IOError:
            fp.close()
            self.Idata = tmp
            raise IOError, 'Could not save the file %s.' % filename

        fp.close()
        self.Idata = tmp

    def len_index(self, l):
        """
	Return the fragment length minus the lower bound of the
	length_range for the fragment. 
	
        Return the length of the B{length_range} not being
	used by the fragment.
	
        @param l: Fragment length.
        """
        return l-self.length_range[0]

    def get_iters(self, l,f):
        """
	Return the list of all iterations for a given pair of fragment 
	length and fragment.
	
	@param l: Fragment length.
	@param f: Current fragment.
	"""
        return self.E[self.len_index(l)][f]

    def _assign_iter(self, l, f, i, hits_dict):
  
        iters = self.get_iters(l,f)
        if len(iters) > i:
            del(iters[i:])
        else:
            i = len(iters)
        iters.append(HitList(hits_dict))

    def _current_iter(self, l, f, i):
        
        iters = self.get_iters(l,f)
        if len(iters) > i:
            del(iters[i:])
        else:
            i = len(iters)
        iters.append(Iteration())
        return iters[i]

    def serial_search(self, jobs, update_func = None):
        """
	Conduct a series of searches specified by jobs.  
	B{update_func} is called after each search.  The 
	number of completed jobs is passed to B{update_func}
	[if provided] after each search.
	
	@param jobs: A dictionary defined in the main GUI. 
	It contains the keys: 'matrix', 'conv_type', 'scale',
	'weights', 'reg', 'iter' and 'cutoff'.  For each key 
	there is a corresponding action (or job) to be completed.
	@param update_func: Number of completed jobs is 
	passed to this function.
	
	"""
	global jobs_counter
        jobs_counter = 0
        I = self.Idata.I
        if update_func != None:
            update_func(jobs_counter)
        for k,v in jobs.items():
            del(jobs[k])
            HL = self._run_search(k, v, I)
            if HL != None:
                self._assign_iter(k[0], k[1], v['iter'], HL)
            jobs_counter += 1
            if update_func != None:
                update_func(jobs_counter)

    def _run_search(self, coords, job, I):
   
        a = coords[1]
        b = a + coords[0]
        qseq = self.rqseq[a:b]
        M, conv_type, PM = _get_converted_matrix(job['matrix'],
                                                job['conv_type'],
                                                job['scale'],
                                                job['weights'],
                                                job['reg'],
                                                self,
                                                list(coords) + [job['iter']])
        r0 = job['cutoff'][1]
        ctype = job['cutoff'][0]

        if ctype == REL_SRCH:
            return None

        if ctype == RNG_SRCH:
            HL=I.rng_srch(qseq, M, r0, conv_type)
        elif ctype == KNN_SRCH:
            HL=I.kNN_srch(qseq, M, r0)

        HL['matrix_name'] = job['matrix']
        HL['matrix'] = PM
        return HL

    def _threaded_search_setup(self, jobs):
        keys = jobs.keys()
        srch_args = []
        for k in keys:
            a = [None]*5
            # query sequence
            a[0] = self.rqseq[k[1]: k[1] + k[0]]

            # matrix - we make a copy for each search
            job = jobs[k]
            M, conv_type, PM = _get_converted_matrix(job['matrix'],
                                                    job['conv_type'],
                                                    job['scale'],
                                                    job['weights'],
                                                    job['reg'],
                                                    self,
                                                    list(k) + [job['iter']])

            jobs[k]['PM'] = PM
            a[1] = M

            # range - 0 if not used
            # kNN - -1 if not used
            ctype = jobs[k]['cutoff'][0]
            if ctype == RNG_SRCH:
                a[2] = jobs[k]['cutoff'][1]
                a[3] = -1
            elif ctype == KNN_SRCH:
                a[2] = 0
                a[3] = jobs[k]['cutoff'][1]
  
            # query defline - leave empty to save memory
            a[4] = ""
            
            # We skip optional arguments - defaults will
            #    be used
            srch_args.append(a)

        return keys, srch_args

    def _threaded_search_assign(self, jobs, keys, res):
        for i,k  in enumerate(keys):
            res[i]['matrix_name'] = jobs[k]['matrix']
            res[i]['matrix'] = jobs[k]['PM']
            self._assign_iter(k[0], k[1], v['iter'], res[i])
            del(jobs[k])

    def threaded_search(self, jobs, I):
        """
	Run the threaded search.  
	
	Each time the search is run a copy of the matrix is made.  
	srch_args is a list which includes the following attributes: 
	the query sequence, the matrix, information on which search 
	was used (0 if range not used, -1 if kNN not used), and the 
	query defline is left blank to save memory.
	
	The search is run recursively.
	"""
        # Construct the list of parameters
        keys, srch_args = self._threaded_search_setup(jobs)
        # Run threaded search
        res = I.threaded_search(srch_args);
        # Get results
        self._threaded_search_assign(jobs, keys, res)

    # Full views for display. Dictionary with seq_id or
    # cluster_id as keys. Values are dictionaries as well
    # with fragments as keys and lists of hits as values.


    def full_view_by_seqid(self, l, iter):
        """
        Full view for display.  
	
	Displayed with fragments in one dimension and hits in the 
	other dimension.
	
	@param l: Fragment length.
	@return: The dictionary and the sorted group.	
	"""
        dict = {}
        for f, frag in enumerate(self.E[self.len_index(l)]):
            if not len(frag): continue 
            if iter >= len(frag): iter = len(frag)-1
            for ht in frag[iter].hits:
                seqid = ht.seq_id
                if seqid not in dict:
                    dict[seqid] = {f:[[ht],0]} 
                elif f not in dict[seqid]:
                    dict[seqid][f] = [[ht],0]
                else:
                    dict[seqid][f][0].append(ht)
                    dict[seqid][f][1] = 1

        # Sort keys by (cluster, seqid)
        grp = [(seqid, self.Idata.seq2cluster(seqid)) for seqid in dict.keys()]
        grp.sort()
        return grp, dict
    
    def full_view_by_cluster(self, l, iter):
        # TO DO: call full_view_by_seqid and then
        # group by cluster
        return self.full_view_by_seqid(l, iter)
    
