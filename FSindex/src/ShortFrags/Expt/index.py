"""
Provides the code to work with index fragment databases.

Only one instance of B{IndexedDb} should be instantiated, 
the B{Idata} (index data) will be set by the most recent index.  

Idata from hit_list.py is assigned dynamically in this module.  
Idata includes the following attributes:
    - B{I}
        - The index.
	- Includes the index name, alphabet, and partitions.
    - B{get_frag}
        - The sequence fragment.
    - B{deflines}
        - Description of the hit.
    - B{accessions}
        - Accession numbers.
    - B{sdb}
        - Pointer to FASTA database.
    
This module may be used to load an index and run a search.  To 
look at the results of a previous search only a FASTA database need 
be loaded.  Primary use of this module is to define index information
B{Idata} for hit_list.

B{Exceptions}:
    - RuntimeError
    - IOError

B{Classes}:
   - IndexedDb  
         -  Class provides the functions to load the databases and 
         print the index.

"""

from cStringIO import StringIO
import os
import os.path
import struct
import md5
import FS
import hit_list
from db import db
from hit_list import HitList

FS_BINS = FS.FS_BINS
SUFFIX_ARRAY = FS.SUFFIX_ARRAY
SEQ_SCAN = FS.SEQ_SCAN

SARRAY = FS.SARRAY
DUPS_ONLY = FS.DUPS_ONLY
FULL_SCAN = FS.FULL_SCAN

def _get_db(I):
    return db(FS.Index_s_db_get(I), new=False, own=False) 

class FSIndex(object):
    def __init__(self, filename, sepn=None, use_sa=1, print_flag=0):
        if sepn == None:
            sepn = []
        self.this = FS.new_Index(filename, sepn, use_sa, print_flag)
        self.thisown = 1
        self.__dict__.update(FS.Index_get_data(self))

    def __del__(self):
        FS.delete_Index(self)

    def save(self, filename):
        return FS.Index_save(self, filename)

    def __str__(self):
        return FS.Index___str__(self)

    def seq2bin(self, seq):
        return FS.Index_seq2bin(self, seq)

    def print_bin(self, bin, options=1):
        return FS.Index_print_bin(self, bin, options)

    def print_stats(self, options=3):
        return FS.Index_print_stats(self, options)

    def get_bin_size(self, bin):
        return FS.Index_get_bin_size(self, bin)

    def get_unique_bin_size(self, bin):
        return FS.Index_get_unique_bin_size(self, bin)

    def rng_srch(self, qseq, M, rng, stype=FS_BINS,
                 ptype=SARRAY, qdef=""):
        
        hits_dict = FS.Index_rng_srch(self, qseq, M, rng, M.conv_type,
                                      stype, ptype, qdef)
        return HitList(hits_dict)

    def kNN_srch(self, qseq, M, kNN, stype=FS_BINS,
                 ptype=SARRAY, qdef=""):
        hits_dict = FS.Index_kNN_srch(self, qseq, M, kNN,
                                      stype, ptype, qdef)
        return HitList(hits_dict)

    def threaded_search(self, srch_args):
        results = FS.Index_threaded_search(self, srch_args)
        return [HitList(r) for r in results]

    db = property(_get_db)


def _md5digest(filename):
 
    BLOCKSIZE = 1024*1024
    fp = open(filename, "rb")
    sum = md5.new()
    while 1:
        block = fp.read(BLOCKSIZE)
        if not block:
            break
        sum.update(block)
    fp.close()
    return sum.hexdigest()


class IndexedDb(object):

    """
    Provides the functions to load the databases and 
    print the index.
    
    B{Exceptions:}
	    - RuntimeError
	    - IOError
    
    """
    def __init__(self):
        """
        Constructor, initializes 
	Idata for hit_list.
        """
        self.I = None
        self.sdb = None
        self.get_frag = lambda i, a, b: 'X'*(b-a)
        
        self.deflines = None
        self.accessions = None

        self._default_deflines()
        hit_list.Hit.Idata = self

    def load_fasta_db(self, filename, md5sum = None):
        """
	Load the FASTA database.
	
	@param filename: Name of the FASTA file.
	@param md5sum:  Compared (if provided) with the md5 digest of 
	the B{filename}.  RuntimeError raised if the two digests do not 
	match. This can help ensure that the database being opened 
	is the same one that was opened previously.
	
	@return: Return the new message digest.
	"""
        # Do _md5digest and compare - throw exception if bad
        new_md5sum = _md5digest(filename)
        if md5sum != None and new_md5sum != md5sum:
            raise RuntimeError, 'Digests for FASTA database do not match.'

        # Clear everything
        if self.deflines == None:
            self._default_deflines()
        self.sdb = None
        self.I = None

        # Load fasta database
        self.sdb = db(filename)
        self.get_frag = self.sdb.get_frag

        # Assign defline attributes
        if self.deflines == None:
            self._fasta_deflines()

        return new_md5sum
        
    def load_index(self, filename, md5sum = None):
        """
	Load the index.
	
	@param filename: Name of the FSindex file.
	@param md5sum:  Compared (if provided) with the md5 digest of 
	the B{filename}.  RuntimeError raised if the two digests do 
	not match. This can help ensure that the database being opened 
	is the same one that was opened previously.
	
	@return: Return the new message digest.
        """
        # Read the name of fasta file - the first 4 bytes of the
        # file contain the length including the trailing '\0'.
        # The string itself follows.
        fp = file(filename, 'rb')
        s = fp.read(4)
        n = struct.unpack('i', s)[0]
        s = fp.read(n-1)
        fp.close()
        # Construct the full path
        head, tail = os.path.split(filename)
        fasta_path = os.path.join(head, s)

        # Do _md5digest and compare - throw exception if bad
        new_md5sum = _md5digest(fasta_path)
        if md5sum != None and new_md5sum != md5sum:
            raise RuntimeError, 'Digests for FASTA database do not match.'

        # Clear everything
        if self.deflines == None:
            self._default_deflines()
        self.sdb = None
        self.I = None

        # Load Index
        try:
            self.I = FSIndex(filename, [])
        except (IOError, RuntimeError):
            raise IOError, 'Could not load index file %s.' % filename
        self.sdb = self.I.db
        self.get_frag = self.sdb.get_frag

        # Assign defline attributes
        if self.deflines == None:
            self._fasta_deflines()

        return new_md5sum

    def load_descriptions(self, filename):
        pass

    def _default_deflines(self):
        self.get_def = lambda i: "Protein Sequence #%d" % i
        self.seq2cluster = lambda i: i
        self.get_accession = lambda i: str(i)

    def _fasta_deflines(self):
        self.get_def = self.sdb.get_def 
        self.seq2cluster = lambda i: i
        self.get_accession = lambda i: str(i)

##     def _clusters_deflines(self):
##         self.get_def = lambda i: self.deflines[i]
##         self.get_accession = lambda i: self.accessions[i]

    def print_str(self):
        """
	Return a printable string with the index's information.
	
	An empty string will be returned if no index is loaded.
	
	@return: Return file string.
	"""
        if self.I == None: return ""

        file_str = StringIO()

        file_str.write("***** Database Details *****\n")
        file_str.write("Database path: %s\n" % self.sdb.db_name)
        file_str.write("Full length: %d\n" % self.sdb.length)
        file_str.write("Number of sequences: %d\n\n" % self.sdb.no_seq)
        file_str.write("***** Index Details *****\n")
        file_str.write("Index path: %s\n" % self.I.index_name)
        file_str.write("Indexed Alphabet: %s\n" % self.I.alphabet)
        file_str.write("Partitions:\n")
        l = len(self.I.ptable)
        for i in range(0,l-2,2):
            file_str.write((("%2.2d. " % i) + self.I.ptable[i] + "    "))
            file_str.write((("%2.2d. " % (i+1)) + self.I.ptable[i+1] + "\n"))

        if l % 2 == 0:
            file_str.write((("%2.2d. " % (l-2)) + self.I.ptable[l-2] + "    "))
        file_str.write((("%2.2d. " % (l-1)) + self.I.ptable[l-1] + "\n"))

        file_str.write("Number of bins: %d\n" % self.I.bins)
        file_str.write("Number of indexed fragments: %d\n" % self.I.fragments)
        file_str.write("Indexed fragment length: %d\n\n" % self.I.indexed_fragment_length)
        
        return file_str.getvalue()

