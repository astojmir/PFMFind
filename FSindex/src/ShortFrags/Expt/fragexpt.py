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
from cStringIO import StringIO
import os, os.path, types, sys, md5, threading, cPickle, gzip

from ShortFrags.Expt.index import IndexedDb
from ShortFrags.Expt.DirichletMix import DirichletMix, freq_counts, henikoff_weights
from ShortFrags.Expt.DirichletMix import BKGRND_PROBS as bg_dict
from ShortFrags.Expt.DirichletInfo import get_mix
from ShortFrags.Expt.matrix import QUASI, MAX, AVG, SCORE, POSITIONAL
from ShortFrags.Expt.matrix import ScoreMatrix, ProfileMatrix
from ShortFrags.Expt.SearchServer import SearchClient


from Bio import SubsMat
from Bio.SubsMat import MatrixInfo
from Bio.SubsMat.FreqTable import FreqTable, FREQ
from Bio.Alphabet import IUPAC
from Bio.Align import Generic, AlignInfo



RNG_SRCH = 0
KNN_SRCH = 1
REL_SRCH = 2

MATRIX_CTYPE = {'None': 0, 'Quasi': QUASI,
                'Avg': AVG, 'Max': MAX}
		

class FullExptDB(object):
    def __init__(self, qseq, qdef='', qrange=None, lrange=None):
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
    
    def get_length(self, l):
        if l >= self.length_range[0]\
           and l <= self.length_range[1]:
            return self.E[l-self.length_range[0]]
        else:
            return None

    def get_frag(self, l, f):
        """
	Return the list of all iterations for a given pair of fragment 
	length and fragment.
	
	@param l: Fragment length.
	@param f: Current fragment.
	"""
        llist = self.get_length(l)
        if llist != None and f >= self.query_range[0]\
               and f < self.query_range[1]:
            return llist[f-self.query_range[0]]
        else:
            return None

    def get_iter(self, l, f, i):
        flist = self.get_frag(l,f)
        if flist != None and i >= 0:
            if i >= len(flist):
                i = -1
            return flist[i]
        else:
            return None
        
    def set_iter(self, l, f, i, HL):
        iters = self.get_frag(l,f)
        if len(iters) > i:
            del(iters[i:])
        else:
            i = len(iters)
        iters.append(HL)
        

class FullExpt(object):
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
        self.Idata = IndexedDb()
        self.search_lock = False
        self.search_client = SearchClient()

        
    def set_params(self, *args):
        self.data = FullExptDB(*args)

    def load_descriptions(self, filename):
        """
        Load the descriptions.  The function initializes B{Idata} by 
        instantiating an instance of B{IndexedDB}. The descriptions
        are then loaded by calling the B{load_descriptions} function
        in B{IndexedDb}.
       
        @param filename: The name of the file containing the descriptions
        to be loaded.
        """
        self.Idata.load_descriptions(filename)


    def load(self, filename):
        """
        Load a saved experiment.
        
        @param filename: The name of the experiment file to load.
        """
        fs = StringIO()
        fp = gzip.GzipFile(filename, 'rb')

        try:
            while 1:
                dt = fp.read()
                if dt == "":
                    break
                fs.write(dt)
            self.data = cPickle.loads(fs.getvalue())
        finally:
            fp.close()
            
    def save(self, filename):
        """
        Save the experiment under the given filename.
        Raise an exception if the file cannot
        be saved.
       
        @param filename: The name the file will be saved under.
        """

        fp = gzip.GzipFile(filename, 'wb')
        try:
            fp.write(cPickle.dumps(self.data, 2))
        finally:
            fp.close()

    def get_matrix(self, matrix_name, conv_type, scale=2.0, weight_type=None,
                   dirichlet_type=None, coords=None):
        PM = None
        matrix_type = SCORE

        if weight_type != None: # Iterated search
            l = coords[0]
            f = coords[1]
            i = coords[2]

            # Assuming that the iteration is correct
            seqs = self.data.get_iter(l,f,i-1).get_seqs()
        
            if conv_type == -1:
                conv_type = 'None'

            # Calculate sequence weights
            bcounts = freq_counts(seqs, "ACDEFGHIKLMNPQRSTVWY")
            if weight_type == 'None':
                weights = [1.0]*len(seqs)
            elif weight_type == 'Henikoff':
                weights = henikoff_weights(seqs, "ACDEFGHIKLMNPQRSTVWY", bcounts)

            if matrix_name == 'PSSM': 
                DM = get_mix(dirichlet_type)
                if weight_type == 'Henikoff':
                    bcounts = DM.block_counts(seqs, weights)
                bprobs = DM.block_probs(bcounts)
                bkgrnd = DM.aa_vector(bg_dict)
                PM = DM.block2pssm(DM.block_log_odds(bprobs, bkgrnd, scale), seqs[0])
                matrix_type = POSITIONAL

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
                PM = lom
        else:
            PM = MatrixInfo.__dict__[matrix_name]

        ctype = MATRIX_CTYPE[conv_type]
        return PM, matrix_type, ctype

    def print_matrix(self, *args):

        PM, matrix_type, ctype = self.get_matrix(*args)
        if matrix_type == SCORE:
            M = ScoreMatrix(PM)
        elif  matrix_type == POSITIONAL:
            M = ProfileMatrix(PM.pssm)
        else:
            return ""

        if ctype:
            M.conv_type = ctype
            M = M.matrix_conv()

        return M.__str__()

    def print_align(self, weight_type, dirichlet_type, coords):
        """
        Return a printable string with the multiple alignment information for 
        the iteration previous to the one given by coords = (l,f,i).
    
        @param weight_type: Type of sequence weights (None or 'Henikoff').
        @param dirichlet_type: Name of Dirichlet mixture.  
        Must be contained in B{NAMES} list in DirichletMix.py.
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
        HL = self.data.get_iter(l,f,i-1)
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
        # file_str.write('\n***** INFORMATION CONTENT *****\n')
        # file_str.write(DM.block2pssm(wcounts, seqs[0]).__str__())
        return file_str.getvalue()

    def search(self, jobs):
        """
	Conduct a series of searches specified by jobs.  
	B{update_func} is called after each search.  The 
	number of completed jobs is passed to B{update_func}
	[if provided] after each search.
	
	@param jobs: A dictionary defined in the main GUI. 
	It contains the keys: 'matrix', 'conv_type', 'scale',
	'weights', 'reg', 'iter' and 'cutoff'.  For each key 
	there is a corresponding action (or job) to be completed.
	"""
        if self.search_lock:
            jobs = {}
            raise RuntimeError,\
                  "Search in progress"
        self.search_lock = True
        
        srch_args = []
        keys = jobs.keys()
        for k in keys:
            job = jobs[k]
            qseq = self.data.rqseq[k[1]: k[1] + k[0]]
            search_type = job['cutoff'][0]
            r0 = job['cutoff'][1]
            PM, matrix_type, conv_type = self.get_matrix(job['matrix'],
                                                         job['conv_type'],
                                                         job['scale'],
                                                         job['weights'],
                                                         job['reg'],
                                                         list(k) + [job['iter']])
            srch_args.append([search_type, qseq, PM,
                                matrix_type, r0, conv_type])

        # Submit search
        results = self.search_client.search(srch_args)
        
        # Process search - assign
        for i,k in enumerate(keys):
            job = jobs[k]
            HL = results[i]
            if HL == None:
                continue
            HL.matrix_name = job['matrix']
            if srch_args[i][3] == POSITIONAL:
                HL.matrix = PM
            else:
                HL.matrix = None
            self.data.set_iter(k[0], k[1], job['iter'], HL)
            del(jobs[k])

        self.search_lock = False

    def full_view_by_seqid(self, l, iter):
        """
        Full view for display.  
	
	Displayed with fragments in one dimension and hits in the 
	other dimension.
	
	@param l: Fragment length.
	@return: The dictionary and the sorted group.	
	"""
        dct = {}
        for f, frag in enumerate(self.data.get_length(l)):
            if not len(frag): continue 
            if iter >= len(frag): iter = len(frag)-1
            for ht in frag[iter]:
                seqid = ht.accession
                if seqid not in dct:
                    dct[seqid] = {f:[[ht],0]} 
                elif f not in dct[seqid]:
                    dct[seqid][f] = [[ht],0]
                else:
                    dct[seqid][f][0].append(ht)
                    dct[seqid][f][1] = 1

        # Sort keys by accession
        grp = [(seqid, seqid) for seqid in dct.keys()]
        grp.sort()
        return grp, dct
    
    def full_view_by_cluster(self, l, iter):
        # TO DO: call full_view_by_seqid and then
        # group by cluster
        return self.full_view_by_seqid(l, iter)

    def keyword_view(self, l, iter):
        dct = {}
        for f, frag in enumerate(self.data.get_length(l)):
            if not len(frag): continue 
            if iter >= len(frag): iter = len(frag)-1
            for ht in frag[iter]:
                for k in ht.keywords:
                    if k not in dct:
                        dct[k] = {f:[[ht],0]} 
                    elif f not in dct[k]:
                        dct[k][f] = [[ht],0]
                    else:
                        dct[k][f][0].append(ht)
                        dct[k][f][1] = 1

        # Sort by keyword (alphabetical)
        grp = [(k,k) for k in dct.keys()]
        grp.sort()
        return grp, dct

