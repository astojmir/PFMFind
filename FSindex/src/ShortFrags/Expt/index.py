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
import os, os.path, struct, re, string, md5
from ShortFrags.Expt import FS, keywords
from ShortFrags.Expt.db import db
from ShortFrags.Expt.hit_list import Hit, HitList



FS_BINS = FS.FS_BINS
SUFFIX_ARRAY = FS.SUFFIX_ARRAY
SEQ_SCAN = FS.SEQ_SCAN

SARRAY = FS.SARRAY
DUPS_ONLY = FS.DUPS_ONLY
FULL_SCAN = FS.FULL_SCAN

def _get_db(I):
    return db(FS.Index_s_db_get(I), new=False, own=False) 

def _get_ix_data(I):
    return FS.Index_get_data(I)

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
    ix_data = property(_get_ix_data)


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
        self.defline = {}
        self.kw = None

    def load_keywords(self, filename):
        self.kw = keywords.SprotKeywords()        
        self.kw.load(filename)
        Hit.keywords = property(fget=lambda hit, kw=self.kw : kw.get_keywords(hit.accession))

    def load_descriptions(self, filename):
        fp = file(filename, 'r')
        for line in fp:
            s = string.split(line)
            self.defline[s[0]] = s[1]
        fp.close()
        Hit.defline = property(fget=lambda hit,\
                               defline=self.defline : defline.get(hit.accession, hit.accession))
 
