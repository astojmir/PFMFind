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

def get_converted_matrix(matrix, conv_type, scale=1.0, weight_type=None,
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
    M, ctype, PM = get_converted_matrix(matrix, conv_type, scale, weight_type,
                                        dirichlet_type, FE, coords)
					
		"""
		Print the converted matrix.
		"""
    return M.__str__()

def print_align(weight_type, dirichlet_type, FE, coords):
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

def filter_ext(s, ext):
    from os.path import splitext
    fsplit = splitext(s)
    return fsplit[1] == ext

def load_FE(filename):
    """
    Load a saved full experiment file.
    
    Exceptions:
        - IOError exception will occur if the file cannot
	be loaded.
    
    @return:  Return the full experiment.
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
    
    Methods:
        - __init__(self)
        - set_params(self, qseq, qdef='', qrange=None, lrange=None)
        - load_fastadb(self, filename)
	- load_index(self, filename)
	- load_descriptions(self, filename)
	- save(self, filename)
	- len_index(self, l)
	- get_iters(self, l,f)
	- assign_iter(self, l, f, i, hits_dict)
	- current_iter(self, l, f, i)
	- serial_search(self, jobs, update_func = None)
	- run_search(self, coords, job, I)
	- threaded_search_setup(self, jobs)
	- threaded_search_assign(self, jobs, keys, res)
	- threaded_search(self, jobs, I)
	- full_view_by_seqid(self, l, iter)
	- full_view_by_cluster(self, l, iter)
    """
    def __init__(self):
        self.mcache = {}
        self.Idata = None
        self.index_filename = None
        self.fasta_filename = None
        self.desc_filename = None
        self.fasta_md5sum = None
        self.query_seq = None
        
    def set_params(self, qseq, qdef='', qrange=None, lrange=None):
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
        if self.Idata == None:
            self.Idata = IndexedDb()
        self.index_filename = None
        self.fasta_filename = None
        self.fasta_md5sum = self.Idata.load_fasta_db(filename,
                                                     self.fasta_md5sum)
        self.fasta_filename = filename
        
    def load_index(self, filename):
        if self.Idata == None:
            self.Idata = IndexedDb()
        self.index_filename = None
        self.fasta_filename = None
        self.fasta_md5sum = self.Idata.load_index(filename,
                                                  self.fasta_md5sum)
        self.index_filename = filename
       
    def load_descriptions(self, filename):
        if self.Idata == None:
            self.Idata = IndexedDb()
        self.desc_filename = None
        self.Idata.load_descriptions(filename)
        self.desc_filename = filename

    def save(self, filename):
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
        return l-self.length_range[0]

    def get_iters(self, l,f):
        return self.E[self.len_index(l)][f]

    def assign_iter(self, l, f, i, hits_dict):
        iters = self.get_iters(l,f)
        if len(iters) > i:
            del(iters[i:])
        else:
            i = len(iters)
        iters.append(HitList(hits_dict))

    def current_iter(self, l, f, i):
        iters = self.get_iters(l,f)
        if len(iters) > i:
            del(iters[i:])
        else:
            i = len(iters)
        iters.append(Iteration())
        return iters[i]

    def serial_search(self, jobs, update_func = None):
        global jobs_counter
        jobs_counter = 0
        I = self.Idata.I
        if update_func != None:
            update_func(jobs_counter)
        for k,v in jobs.items():
            del(jobs[k])
            HL = self.run_search(k, v, I)
            if HL != None:
                self.assign_iter(k[0], k[1], v['iter'], HL)
            jobs_counter += 1
            if update_func != None:
                update_func(jobs_counter)

    def run_search(self, coords, job, I):
        a = coords[1]
        b = a + coords[0]
        qseq = self.rqseq[a:b]
        M, conv_type, PM = get_converted_matrix(job['matrix'],
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

    def threaded_search_setup(self, jobs):
        keys = jobs.keys()
        srch_args = []
        for k in keys:
            a = [None]*5
            # query sequence
            a[0] = self.rqseq[k[1]: k[1] + k[0]]

            # matrix - we make a copy for each search
            job = jobs[k]
            M, conv_type, PM = get_converted_matrix(job['matrix'],
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

    def threaded_search_assign(self, jobs, keys, res):
        for i,k  in enumerate(keys):
            res[i]['matrix_name'] = jobs[k]['matrix']
            res[i]['matrix'] = jobs[k]['PM']
            self.assign_iter(k[0], k[1], v['iter'], res[i])
            del(jobs[k])

    def threaded_search(self, jobs, I):
        # Construct the list of parameters
        keys, srch_args = self.threaded_search_setup(jobs)
        # Run threaded search
        res = I.threaded_search(srch_args);
        # Get results
        self.threaded_search_assign(jobs, keys, res)

    # Full views for display. Dictionary with seq_id or
    # cluster_id as keys. Values are dictionaries as well
    # with fragments as keys and lists of hits as values.


    def full_view_by_seqid(self, l, iter):
        """
	Display 
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
    
