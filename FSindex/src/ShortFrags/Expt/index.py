"""
This module loads and prints the index.

Exceptions:
    - RuntimeError
    - IOError

Classes:
   - IndexedDb  
         -  Load the databases and print the index.

Functions:
   - md5digest(filename) 
         - Read a file in as a message digest.


"""

from cStringIO import StringIO
import os
import os.path
import struct
import md5
from ShortFrags import FS
import hit_list


BLOCKSIZE = 1024*1024
def md5digest(filename):
    """
    Open a file and read it in as a new message digest.

    @return: Return a message digest as a string of length 
    32 with only hexadecimal digits.
    """
    fp = open(filename, "rb")
    sum = md5.new()
    while 1:
        block = fp.read(BLOCKSIZE)
        if not block:
            break
        sum.update(block)
    fp.close()
    return sum.hexdigest()


class IndexedDb:

    """
    Class contains the functions to load the databases and 
    print the index.
    
    Methods:
        - __init__(self)
        - load_fasta_db(self, filename, md5sum = None)
        - load_index(self, filename, md5sum = None)
        - load_descriptions(self, filename)
        - print_str(self)
    """
    def __init__(self):
        self.I = None
        self.sdb = None
        self.get_frag = lambda i, a, b: 'X'*(b-a)
        
        self.deflines = None
        self.seq2clusters = None
        self.clusters2seqs = None
        self.accessions = None

        self._default_deflines()
        hit_list.Hit.Idata = self

    def load_fasta_db(self, filename, md5sum = None):
        """
	Load the fasta database.
	
	Exceptions:
	    - Throws a RuntimeError exception if the loaded 
	    database does not match the database given as 
	    an argument.
	
	@return: Return the new message digest.
	"""
        # Do md5digest and compare - throw exception if bad
        new_md5sum = md5digest(filename)
        if md5sum != None and new_md5sum != md5sum:
            raise RuntimeError, 'Digests for FASTA database do not match.'

        # Clear everything
        if self.deflines == None:
            self._default_deflines()
        self.sdb = None
        self.I = None

        # Load fasta database
        self.sdb = FS.db(filename)
        self.get_frag = self.sdb.get_frag

        # Assign defline attributes
        if self.deflines == None:
            self._fasta_deflines()

        return new_md5sum
        
    def load_index(self, filename, md5sum = None):
        """
	Load the index.
	
	Exceptions: 
	    - Throws a RuntimeError exception if the loaded 
	    database does not match the fasta database.
	    - Throws an IOError if the index file cannot be loaded.
	
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

        # Do md5digest and compare - throw exception if bad
        new_md5sum = md5digest(fasta_path)
        if md5sum != None and new_md5sum != md5sum:
            raise RuntimeError, 'Digests for FASTA database do not match.'

        # Clear everything
        if self.deflines == None:
            self._default_deflines()
        self.sdb = None
        self.I = None

        # Load Index
        try:
            self.I = FS.index(filename, [])
        except (IOError, RuntimeError):
            raise IOError, 'Could not load index file %s.' % filename
        self.sdb = self.I.s_db
        self.get_frag = self.sdb.get_frag

        # Assign defline attributes
        if self.deflines == None:
            self._fasta_deflines()

        return new_md5sum

    def load_descriptions(self, filename):
        """
	Does nothing, method employs pass statement.
	"""
        pass

    def _default_deflines(self):
        self.get_def = lambda i: "Protein Sequence #%d" % i
        self.seq2cluster = lambda i: i
        self.get_accession = lambda i: str(i)

    def _fasta_deflines(self):
        self.get_def = self.sdb.get_def 
        self.seq2cluster = lambda i: i
        self.get_accession = lambda i: str(i)

    def _clusters_deflines(self):
        self.get_def = lambda i: self.deflines[i]
        self.seq2cluster = lambda i: self.seq2clusters[i]
        self.get_accession = lambda i: self.accessions[i]

    def print_str(self):
        """
	Print the index's information.
	
	Will print nothing if index does not exist.
	
	@return: Return file string.
	"""
        if self.I == None: return ""

        class Data:
            pass
        data = Data()
        data.__dict__ = self.I.get_data()
        file_str = StringIO()

        file_str.write("***** Database Details *****\n")
        file_str.write("Database path: %s\n" % data.db_name)
        file_str.write("Full length: %d\n" % data.db_length)
        file_str.write("Number of sequences: %d\n\n" % data.db_no_seq)
        file_str.write("***** Index Details *****\n")
        file_str.write("Index path: %s\n" % data.index_name)
        file_str.write("Indexed Alphabet: %s\n" % data.alphabet)
        file_str.write("Partitions:\n")
        l = len(data.ptable)
        for i in range(0,l-2,2):
            file_str.write((("%2.2d. " % i) + data.ptable[i] + "    "))
            file_str.write((("%2.2d. " % (i+1)) + data.ptable[i+1] + "\n"))

        if l % 2 == 0:
            file_str.write((("%2.2d. " % (l-2)) + data.ptable[l-2] + "    "))
        file_str.write((("%2.2d. " % (l-1)) + data.ptable[l-1] + "\n"))

        file_str.write("Number of bins: %d\n" % data.bins)
        file_str.write("Number of indexed fragments: %d\n" % data.fragments)
        file_str.write("Indexed fragment length: %d\n\n" % data.indexed_fragment_length)
        
        return file_str.getvalue()

